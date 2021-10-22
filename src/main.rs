extern crate csv;
extern crate ndarray_csv;

use bio::io::bed;
use structopt::StructOpt;
use std::convert::TryInto;
//use ndarray::prelude::*;
use ndarray_csv::Array2Writer;
use rust_htslib::{bam, bam::Read};

use vplot::{VMatrix, VEntry, VRegions};

#[derive(StructOpt, Debug)]
#[structopt(setting = structopt::clap::AppSettings::ArgRequiredElseHelp)]
struct Cli {
    /// Path to an indexed bam file
    #[structopt(parse(from_os_str))]
    bam: std::path::PathBuf,
    /// Path to a bed file (must be in bed4 format: chr, start, end, strand) Of
    /// a region (or regions) in which to generate the vplot. If using multiple
    /// regions, all entries must be the same width. If setting multiple
    /// regions, reads will be aggregated into a single matrix unless `--multi`
    /// is set.
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

    let mut ventry = VEntry::new();

    let mut vmatrix = VMatrix::new(&regions, &args.max_fragment_size, &args.multi, &args.invert);

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
                    "midpoint" => vmatrix.insert_midpoint(&ventry),
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

    fn get_output_filename(path: &str, bam: &str) -> String {
        // if path is directory write to dir/vmatrix.csv
        // https://stackoverflow.com/a/64952550

        let outfile;

        if path.ends_with("/") {
            outfile = format!("{}/{}.vmatrix.csv", &path, &bam);
        } else {
            outfile = path.to_string();
        }

        outfile

    }

    let output_path : String = get_output_filename(output_prefix, &args.bam.into_os_string().into_string().unwrap());

    if vmatrix.aggregate {
        match output_prefix {
            "" => vmatrix.write_csv_stdout(),
             _ => vmatrix.write_csv(&output_path)
        }
    }
}
