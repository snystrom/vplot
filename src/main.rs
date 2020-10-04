extern crate bio;
extern crate matrix;

use std::convert::{TryInto, TryFrom};
use rust_htslib::{bam, bam::Read};
//use rust_htslib::bam;
use structopt::StructOpt;
use bio::io::bed;
use matrix::prelude::*;

#[derive(StructOpt, Debug)]
#[structopt(setting = structopt::clap::AppSettings::ArgRequiredElseHelp)]
struct Cli {
    #[structopt(parse(from_os_str))]
    bam: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    regions: std::path::PathBuf,
    /// Maximum fragment size to include in the V-plot matrix
    #[structopt(default_value = "700", short = "x", long = "max-size")]
    max_fragment_size: i64,

}

struct VEntry {
    /// Start position of fragment
    start: i64,
    /// Insert size of the fragment
    insert_size: i64,
    /// End coordinate of the bed region
    region_end: i64,
    /// index of which bed region this fragment resides in
    // must be i64 because it's used as a multiplier for other i64s
    // so, future me, don't try getting clever
    region_n: i64,
}

impl VEntry {
    /// Initialize struct with -1's
    fn allocate() -> VEntry {

        VEntry {
            start: -1,
            insert_size: -1,
            region_end: -1,
            region_n: -1,
        }
    }

    /// Midpoint of the fragment
    fn midpoint(&self) -> i64 {
        // NOTE: midpoint is rounded down if insert_size is odd
        // TODO: Fix?? Round up?
        self.start + (self.insert_size/2)
    }

    /// Row in the VMatrix this entry should be inserted into
    fn row(&self) -> usize {
        // NOTE: region count has to be > 0, else will get all 0's for first region row positions
        ((self.region_n + 1) * self.insert_size.abs()).try_into().unwrap()
    }

    /// Column in the VMatrix this entry should be inserted into
    // value must be usize because matrix library needs usize
    fn col(&self) -> usize {
        (self.region_end - self.midpoint()).try_into().unwrap()
    }

    /// Update fragment entry
    /// Allows operation in same allocation
    fn update(&mut self, region: &bed::Record, region_n: i64, record: &bam::Record) {
        self.start = record.pos();
        self.insert_size = record.insert_size();
        self.region_end = region.end().try_into().unwrap();
        self.region_n = region_n;
    }
}

// TODO: VMatrix struct that hold matrix

fn count_bed_entries(path: &std::path::PathBuf) -> i64 {
    let bed_reader = bed::Reader::from_file(path);
    bed_reader.unwrap()
              .records()
              .count() as i64
}

/// Return (n_bed_entries, bed_width)
/// n_bed_entries * max_fragment_size = nrow
fn get_dims_from_bed(path: &std::path::PathBuf) -> (usize, usize) {
    // loop over records and ensure same width
    // count records

    let bed_reader = bed::Reader::from_file(path);
    let mut n_regions = 0;
    let mut set_width = 0;
    let mut init = true;
    for region in bed_reader.unwrap().records() {
        let record = region.unwrap();
        n_regions += 1;
        // Grab width of first region
        // all other regions should be this width
        if init {
            set_width = record.end() - record.start();
            init = false;
        }

        let width = record.end() - record.start();

        if set_width != width {
            panic!("All bed regions must be equal width\nWidth: {} of bed entry on line {} does not match width of first entry: {}", width, n_regions, set_width);
        }

    }

    (n_regions as usize, set_width as usize)
}

/// (nrow, ncol)
fn get_matrix_dims(args: &Cli) -> (usize, usize){
    let (n_bed, ncol) = get_dims_from_bed(&args.regions);
    let mfs = usize::try_from(args.max_fragment_size).unwrap();
    let nrow = n_bed * mfs;
    (nrow, ncol)
}

fn connect_indexed_bam(path: &std::path::PathBuf) -> bam::IndexedReader {
    bam::IndexedReader::from_path(path)
        .ok()
        .expect("Cannot read bam file")
}

fn initialize_matrix(nrow: usize, ncol: usize) -> Compressed<usize> {
    Compressed::zero((nrow, ncol))
}

fn main() {
    let args = Cli::from_args();
    let bed_reader = bed::Reader::from_file(&args.regions);
    let mut bam_reader = connect_indexed_bam(&args.bam);

    let mut ventry = VEntry::allocate();

    // TODO: remove below comments (debug)
    //let (n_bed, ncol) = get_dims_from_bed(&args.regions);
    //let nrow = n_bed * (args.max_fragment_size as usize);
    let (nrow, ncol) = get_matrix_dims(&args);

    let mut vmatrix = initialize_matrix(nrow, ncol);

    let mut bed_region_n = 0 as i64;
    for bed_region in bed_reader.unwrap().records() {
        let region = bed_region.ok()
                               .expect("Error reading bed region");

        //println!("{}\t{}\t{}", region.chrom(), region.start(), region.end());

        let region_chrom = region.chrom();
        let bed_range = std::ops::Range{start: region.start(), end: region.end()};
        let tid = bam_reader.header().tid(region_chrom.as_bytes()).unwrap();

        bam_reader.fetch(tid, region.start(), region.end()).unwrap();

        // read fetched records into bam_record in loop,
        // write out to matrix
        let mut record = bam::Record::new();

        // initialize loop
        // check holds bool determining if reads remain in region
        let mut check = bam_reader.read(&mut record);

        //TODO:
        // Since I want an indexed bam file, it can't be sorted.
        // so, I think what I need to do is for each record
        // only save if record.is_paired(), on the same chromosome, and totally within the bed region, and insert size < size cutoff.
        //  -
        // Then save each midpoint of the fragment as: record.pos() + record.insert()
        while check.unwrap() {

            // TODO: if record is totally within region && insert < cutoff
            // NOTE: Actually, I think only midpoint needs to be within the region, not the whole read.
            if record.is_proper_pair() && record.is_first_in_template() {
                ventry.update(&region, bed_region_n, &record);
                if bed_range.contains(&ventry.midpoint().try_into().unwrap()) && ventry.insert_size.abs() <= args.max_fragment_size {
                    // TODO: write to matrix
                    let mut val = vmatrix.get((ventry.row(), ventry.col()));
                    val += 1;
                    vmatrix.set((ventry.row(), ventry.col()), val);
                    //println!("{}\t{}\t{}", ventry.row(), ventry.col(), ventry.midpoint());
                }
            }

            // Advance to next bam region
            check = bam_reader.read(&mut record);
        }

        // Advance to next bed region
        bed_region_n += 1;
    }

    // TODO: remove, print final matrix
    println!("{:?}", vmatrix);
}
