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
    fragment_type: String,
    /// Invert the matrix so that the smallest fragments appear at the top.
    #[structopt(short = "i", long)]
    invert: bool,
    /// Instead of aggregating reads into 1 matrix, write 1 matrix for each region.
    /// Matrices are written as 1 csv per region named: `chr-start-end.csv`
    #[structopt(short = "m", long)]
    multi: bool,
    /// Set output file name or output directory. This option behaves differently depending on which input flags are set. See --help for details.
    ///
    /// If --multi is unset and -o is set to a directory, the output file will be written to: outdir/<bamfile>.vmatrix.csv.
    /// if --multi is unset and -o is a file path, output file will be written to this file name.
    /// if --multi is set and -o is a directory, files will be written to outdir as: outdir/chr-start-end.csv.
    /// if --multi is set and -o is a string, the string will be used as a prefix, and files will be written as: <prefix>chr-start-end.csv.
    ///
    /// Examples:
    ///
    /// vplot reads.bam regions.bed > output.csv
    ///
    /// vplot -o outdir/ reads.bam regions.bed
    ///
    ///   returns: outdir/reads.bam.vmatrix.csv
    ///
    /// vplot -o matrix.csv reads.bam regions.bed
    ///
    ///   returns: matrix.csv
    ///
    /// vplot -m -o outdir/ reads.bam regions.bed
    ///
    ///   returns:
    ///      - outdir/chr1-1000-2000.csv
    ///     
    ///      - outdir/chr2-1000-2000.csv
    ///
    /// vplot -m -o myPrefix_ reads.bam regions.bed
    ///
    ///   returns:
    ///
    ///      - myPrefix_chr1-1000-2000.csv
    ///
    ///      - myPrefix_chr2-1000-2000.csv
    ///
    /// vplot -m -o outdir/myPrefix_ reads.bam regions.bed
    ///
    ///   returns:
    ///
    ///      - outdir/myPrefix_chr1-1000-2000.csv
    ///
    ///      - outdir/myPrefix_chr2-1000-2000.csv
    ///
    #[structopt(default_value = "-", short = "o", long)]
    output: String,
    //TODO: add flag to use bed column as matrix output name

}

struct VEntry {
    /// Start position of fragment
    start: i64,
    /// Insert size of the fragment
    insert_size: i64,
    /// Start coordinate of the bed region
    region_start: i64,
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
            /// Start point of fragment
            region_start: -1,
            /// End point of fragment
            region_end: -1,
            /// Tracks which bed region the read belongs to
            region_n: -1
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
        self.region_start = region.start().try_into().unwrap();
        self.region_end = region.end().try_into().unwrap();
        self.region_n = region_n

    }

    fn row(&self) -> usize {
        // Row in the VMatrix this entry should be inserted into
        // NOTE: region count has to be > 0, else will get all 0's for first region row positions
        // But ndarray is 0 indexed, so subtract 1 from final to get row #
        //((self.insert_size.abs()) - 1).try_into().unwrap()
        (self.insert_size.abs()).try_into().unwrap()
    }
}

/// Holds the vplot matrix
struct VMatrix {
    matrix: ndarray::Array2<i64>,
    ncol: i64,
    //nrow: i64,
    max_fragment_size: i64,
    aggregate: bool,
    /// Whether entries should be inserted using inverted y-coordinates
    invert: bool
}

impl VMatrix {
    fn allocate(regions: &VRegions, max_fragment_size: &i64, multi: &bool, invert: &bool) -> VMatrix {

        let nrow = *max_fragment_size;

        VMatrix {
            matrix: VMatrix::_initialize_matrix(nrow, regions.width),
            ncol: regions.width,
            //nrow: nrow,
            max_fragment_size: *max_fragment_size,
            aggregate: !multi,
            invert: *invert
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
            let col = entry.midpoint() - entry.region_start;

            let row = if self.invert {entry.row()} else {self.max_fragment_size as usize - entry.row()};

            //self.matrix[[entry.row(), col as usize]] += 1;
            self.matrix[[row, col as usize]] += 1;
        }
    }

    fn insert_fragment_ends(&mut self, entry: &VEntry) {
        if entry.insert_size.abs() <= self.max_fragment_size {

            // since I don't care about directionality,
            // i can swap start/end coords of the read fine.
            let start_col = entry.start - entry.region_start;
            let end_col = (entry.start + entry.insert_size) - entry.region_start;

            let row = if self.invert {entry.row()} else {self.max_fragment_size as usize - entry.row()};
            // Possible valid column values
            let col_range = 0..self.ncol;

            // since orientation doesn't matter, just position,
            // use min/max to convert into following tuple:
            // (start, end)
            let (start, end) = (core::cmp::min(start_col, end_col), core::cmp::max(start_col, end_col));

            if col_range.contains(&start) {
                self.matrix[[row, start as usize]] += 1;
            }

            if col_range.contains(&end) {
                self.matrix[[row, end as usize]] += 1;
            }

        }
    }

    fn insert_fragment(&mut self, entry: &VEntry) {
        if entry.insert_size.abs() <= self.max_fragment_size {

            //let start_col = entry.region_end - entry.start;
            //let end_col = entry.region_end - (entry.start + entry.insert_size);
            let start_col = entry.start - entry.region_start;
            let end_col = (entry.start + entry.insert_size) - entry.region_start;

            let row = if self.invert {entry.row()} else {self.max_fragment_size as usize - entry.row()};

            // Possible valid column values
            let col_range = 0..self.ncol;
            // column coordinates of fragment bases
            let fragment_range = core::cmp::min(start_col, end_col)..core::cmp::max(start_col, end_col);

            for col in fragment_range {
                if col_range.contains(&col) {
                    self.matrix[[row, col as usize]] += 1;
                }
            }

        }

    }

    // re-zeroes matrix
    fn clear(&mut self) {
        self.matrix.iter_mut().for_each(|m| *m = 0)
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
            path,
            n_regions,
            width
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
            //eprintln!("{:?}", region);
            let record = region.unwrap();
            n_regions += 1;
            // Grab width of first region
            // all other regions should be this width
            if init {
                //eprintln!("{:?}", record);
                // TODO: throw error if set.width is negative. should end always be > start?
                // maybe also consider parsing the strand column & conditionally enforcing this. ie if "." take abs, if + do end - start, if - do start - end
                set_width = (record.end() - record.start()).try_into().unwrap();
                init = false;
            }

            //TODO: throw error if set.width is negative. should end always be > start?
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

fn get_output_prefix(path: &str) -> &str {
    match path {
        "-" => "",
        _ => path
    }
}

fn main() {
    let args = Cli::from_args();
    let regions = VRegions::new(&args.regions);
    let bed_reader = bed::Reader::from_file(regions.path);
    let mut bam_reader = connect_indexed_bam(&args.bam);
    let output_prefix = get_output_prefix(&args.output);

    let mut ventry = VEntry::allocate();

    let mut vmatrix = VMatrix::allocate(&regions, &args.max_fragment_size, &args.multi, &args.invert);

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
        // TODO: if -multi = True, increment this, else skip (to aggregate into 1 matrix)
        if !vmatrix.aggregate {
            bed_region_n += 1;

            // write csv named "{chr}-{start}-{end}.csv"
            // then re-zero matrix
            let region_name = format!("{}{}-{}-{}.csv", output_prefix, region_chrom, region.start(), region.end());

            let mut writer = csv::Writer::from_path(&region_name)
                .expect(&format!("cannot open connection to {} for writing", &region_name));

            writer.serialize_array2(&vmatrix.matrix)
                  .expect(&format!("cannot write matrix to {}", &region_name));

            vmatrix.clear();

        }
    }


    fn write_mat_stdout(vmatrix: &VMatrix){
        let mut writer = csv::Writer::from_writer(std::io::stdout());
        writer.serialize_array2(&vmatrix.matrix)
              .expect("cannot write matrix to stdout");

    }

    fn write_mat_path(vmatrix: &VMatrix, path: &str, bam: &str){

        // if path is directory write to dir/vmatrix.csv
        // https://stackoverflow.com/a/64952550

        let outfile;

        if path.ends_with("/") {
            outfile = format!("{}/{}.vmatrix.csv", &path, &bam);
        } else {
            outfile = path.to_string();
        }

        let mut writer = csv::Writer::from_path(&outfile)
            .expect(&format!("cannot open connection to {} for writing", &outfile));

        writer.serialize_array2(&vmatrix.matrix)
              .expect(&format!("cannot write matrix to {}", &outfile));

    }

    if vmatrix.aggregate {
        match output_prefix {
            "" => write_mat_stdout(&vmatrix),
            _  => write_mat_path(&vmatrix, output_prefix, &args.bam.into_os_string().into_string().unwrap())
        }
    }
}
