extern crate bio;
extern crate csv;
extern crate ndarray;
extern crate ndarray_csv;

use std::convert::TryInto;
use rust_htslib::{bam, bam::Read};
use structopt::StructOpt;
use bio::io::bed;
use ndarray::prelude::*;
use ndarray_csv::Array2Writer;

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
            /// Start coordinate of read
            start: -1,
            /// Total fragment length (TLEN)
            insert_size: -1,
            /// End point of fragment
            region_end: -1,
            /// Tracks which bed region the read belongs to
            region_n: -1,
        }
    }

    /// Computes the midpoint of the fragment
    fn midpoint(&self) -> i64 {
        // NOTE: midpoint is rounded down if insert_size is odd
        // TODO: Fix?? Round up?
        self.start + (self.insert_size/2)
    }

    /// Row in the VMatrix this entry should be inserted into
    fn row(&self) -> usize {
        // NOTE: region count has to be > 0, else will get all 0's for first region row positions
        // But ndarray is 0 indexed, so subtract 1 from final to get row #
        (((self.region_n + 1) * self.insert_size.abs()) - 1).try_into().unwrap()
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

/// Holds the vplot matrix
struct VMatrix {
    matrix: ndarray::Array2<i64>,
    n_regions: i64,
    width: i64,
    max_fragment_size: i64
}

impl VMatrix {
    fn allocate(n_regions: i64, width: i64, max_fragment_size: i64) -> VMatrix {
        VMatrix{
            matrix: initialize_matrix(width, n_regions),
            n_regions: n_regions,
            max_fragment_size: max_fragment_size
        }


    }
    fn insert_midpoint() {

    }
}

struct VRegions<'a> {
    path: &'a std::path::PathBuf,
    n_regions: i64,
    width: i64,
}

impl VRegions<'_> {
    fn new(path: &std::path::PathBuf) -> VRegions {
        let (n_regions, width) = VRegions::_get_dims_from_bed(path);

        VRegions {
            path: path,
            n_regions: n_regions,
            width: width
        }
    }

    /// Return (n_bed_entries, bed_width)
    /// n_bed_entries * max_fragment_size = nrow
    fn _get_dims_from_bed(path: &std::path::PathBuf) -> (i64, i64) {
        // loop over records and ensure same width
        // count records

        let bed_reader = bed::Reader::from_file(path);
        let mut n_regions = 0 as i64;
        let mut set_width = 0 as i64;
        let mut init = true;
        for region in bed_reader.unwrap().records() {
            let record = region.unwrap();
            n_regions += 1;
            // Grab width of first region
            // all other regions should be this width
            if init {
                set_width = (record.end() - record.start()).try_into().unwrap();
                init = false;
            }

            let width = (record.end() - record.start()).try_into().unwrap();

            if set_width != width {
                panic!("All bed regions must be equal width\nWidth: {} of bed entry on line {} does not match width of first entry: {}", width, n_regions, set_width);
            }

        }
        (n_regions, set_width)
    }

}


/// (nrow, ncol)
fn get_matrix_dims(args: &Cli) -> (i64, i64){
    let (n_bed, ncol) = get_dims_from_bed(&args.regions);
    let nrow = n_bed * args.max_fragment_size;
    (nrow, ncol)
}

fn connect_indexed_bam(path: &std::path::PathBuf) -> bam::IndexedReader {
    bam::IndexedReader::from_path(path)
        .ok()
        .expect("Cannot read bam file")
}

fn initialize_matrix(nrow: i64, ncol: i64) -> Array2<i64> {
    Array2::<i64>::zeros((nrow as usize, ncol as usize))
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

            // NOTE: Here I'm parsing only 0x040 which corresponds to R1 in
            // illumina sequencing This lets me ignore the problem of
            // double-counting read pairs, because the information I need is
            // stored in R1 and R2, I ignore R2.
            if record.is_proper_pair() && record.is_first_in_template() {
                ventry.update(&region, bed_region_n, &record);
                // If midpoint is within the region and below cutoff, save it
                if bed_range.contains(&ventry.midpoint().try_into().unwrap()) && ventry.insert_size.abs() <= args.max_fragment_size {
                    vmatrix[[ventry.row(), ventry.col()]] += 1;
                }
            }

            // Advance to next bam region
            check = bam_reader.read(&mut record);
        }

        // Advance to next bed region
        bed_region_n += 1;
    }

    // TODO: remove, print final matrix
    //println!("{:?}", vmatrix);

    //let file = File::create("mat.csv");
    //let mut writer = WriterBuilder::new() .has_headers(false) .from_writer(std::io::stdout());
    let mut writer = csv::Writer::from_writer(std::io::stdout());
    writer.serialize_array2(&vmatrix)
          .expect("cannot write matrix to stdout");

}
