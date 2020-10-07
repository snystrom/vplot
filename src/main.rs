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
    /// Path to an indexed bam file
    #[structopt(parse(from_os_str))]
    bam: std::path::PathBuf,
    /// Path to a bed file (must be in bed4 format: chr, start, end, strand) Of
    /// a region (or regions) in which to generate the vplot. If using multiple
    /// regions, all entries must be the same width. (NOTE: using more than 1 region is poorly supported currently, but technically works)
    #[structopt(parse(from_os_str))]
    regions: std::path::PathBuf,
    /// Maximum fragment size to include in the V-plot matrix
    #[structopt(default_value = "700", short = "x", long = "max-size")]
    max_fragment_size: i64,
    /// How reads are counted in the matrix. Using either the midpoint of the fragment, fragment ends, or the whole fragment.
    #[structopt(default_value = "midpoint", short = "f", long = "fragment-type", possible_values = &["midpoint", "ends", "fragment"])]
    fragment_type: String

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

    /// Update fragment entry
    /// Allows operation in same allocation
    fn update(&mut self, region: &bed::Record, region_n: i64, read: &bam::Record) {
        self.start = read.pos();
        self.insert_size = read.insert_size();
        self.region_end = region.end().try_into().unwrap();
        self.region_n = region_n

    }

    fn row(&self) -> usize {
        // Row in the VMatrix this entry should be inserted into
        // NOTE: region count has to be > 0, else will get all 0's for first region row positions
        // But ndarray is 0 indexed, so subtract 1 from final to get row #
        (((self.region_n + 1) * self.insert_size.abs()) - 1).try_into().unwrap()
    }
}

/// Holds the vplot matrix
struct VMatrix {
    matrix: ndarray::Array2<i64>,
    n_regions: i64,
    ncol: i64,
    nrow: i64,
    max_fragment_size: i64
}

impl VMatrix {
    fn allocate(regions: &VRegions, max_fragment_size: &i64) -> VMatrix {

        let nrow: i64 = regions.n_regions * max_fragment_size;

        VMatrix {
            matrix: VMatrix::_initialize_matrix(nrow, regions.width),
            ncol: regions.width,
            nrow: nrow,
            n_regions: regions.n_regions,
            max_fragment_size: *max_fragment_size
        }

    }

    fn _initialize_matrix(nrow: i64, ncol: i64) -> Array2<i64> {
        Array2::<i64>::zeros((nrow as usize, ncol as usize))
    }

    fn insert_midpoint(&mut self, entry: &VEntry) {
        //TODO: consider doing the ncol check as in insert_fragment_ends?
        if entry.insert_size.abs() <= self.max_fragment_size {
            // Column in the VMatrix this entry should be inserted into
            // value must be usize because matrix library needs usize
            let col: usize = (entry.region_end - entry.midpoint()).try_into().unwrap();

            self.matrix[[entry.row(), col]] += 1;
        }
    }

    fn insert_fragment_ends(&mut self, entry: &VEntry) {
        if entry.insert_size.abs() <= self.max_fragment_size {

            let start_col = entry.region_end - entry.start;
            let end_col = entry.region_end - (entry.start + entry.insert_size);

            if start_col >= 0 && start_col <= self.ncol as i64 {
                self.matrix[[entry.row(), start_col as usize]] += 1;
            }

            if end_col >= 0 && end_col <= self.ncol as i64 {
                self.matrix[[entry.row(), end_col as usize]] += 1;
            }
        }
    }

    fn insert_fragment(&mut self, entry: &VEntry) {
        if entry.insert_size.abs() <= self.max_fragment_size {

            let start_col = entry.region_end - entry.start;
            let end_col = entry.region_end - (entry.start + entry.insert_size);

            // Possible valid column values
            let col_range = 0..self.ncol;
            // column coordinates of fragment bases
            let fragment_range = core::cmp::min(start_col, end_col)..core::cmp::max(start_col, end_col);

            for col in fragment_range {
                if col_range.contains(&col) {
                    self.matrix[[entry.row(), col as usize]] += 1;
                }
            }

        }

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


fn connect_indexed_bam(path: &std::path::PathBuf) -> bam::IndexedReader {
    bam::IndexedReader::from_path(path)
        .ok()
        .expect("Cannot read bam file")
}

fn main() {
    let args = Cli::from_args();
    let regions = VRegions::new(&args.regions);
    let bed_reader = bed::Reader::from_file(regions.path);
    let mut bam_reader = connect_indexed_bam(&args.bam);

    let mut ventry = VEntry::allocate();

    // TODO: remove below comments (debug)
    //let (n_bed, ncol) = get_dims_from_bed(&args.regions);
    //let nrow = n_bed * (args.max_fragment_size as usize);
    // NOTE: OLD
    //let (nrow, ncol) = get_matrix_dims(&args);
    //let mut vmatrix = initialize_matrix(nrow, ncol);
    let mut vmatrix = VMatrix::allocate(&regions, &args.max_fragment_size);

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
        let mut read = bam::Record::new();

        // initialize loop
        // check holds bool determining if reads remain in region
        let mut check = bam_reader.read(&mut read);

        while check.unwrap() {

            // NOTE: Here I'm parsing only 0x040 which corresponds to R1 in
            // illumina sequencing This lets me ignore the problem of
            // double-counting read pairs, because the information I need is
            // stored in R1 and R2, I ignore R2.
            if read.is_proper_pair() && read.is_first_in_template() {
                ventry.update(&region, bed_region_n, &read);

                match &args.fragment_type[..] {
                    "ends" => vmatrix.insert_fragment_ends(&ventry),
                    "fragment" => vmatrix.insert_fragment(&ventry),
                    "midpoint" => {
                        // TODO: move this logic into the insert function like the other two.
                        if bed_range.contains(&ventry.midpoint().try_into().unwrap()) {
                            vmatrix.insert_midpoint(&ventry);
                        }
                    },
                    _ => println!("Error: Unknown fragment type"),

                    }

                } 

            // Advance to next bam region
            check = bam_reader.read(&mut read);

            }

        // Advance to next bed region
        bed_region_n += 1;
    }

    let mut writer = csv::Writer::from_writer(std::io::stdout());
    writer.serialize_array2(&vmatrix.matrix)
          .expect("cannot write matrix to stdout");

}
