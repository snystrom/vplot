extern crate csv;
extern crate ndarray_csv;

use bio::io::bed;
use structopt::StructOpt;
use rust_htslib::{bam, bam::Read};

use vplot::{VMatrix, VEntry, VRegions};

use plotly::common::{ColorScale, ColorScalePalette, Title};
use plotly::contour::Contours;
use plotly::{Contour, HeatMap, Layout, Plot};

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
    /// Matrices are written as 1 csv per region named: `chr-start-end.vmatrix.csv`
    #[structopt(short = "m", long)]
    multi: bool,
    /// Set output file name or output directory. This option behaves differently depending on which input flags are set. See --help for details.
    ///
    /// If --multi is unset and -o is set to a directory, the output file will be written to: outdir/<bamfile>.vmatrix.csv.
    /// if --multi is unset and -o is a file path, output file will be written to this file name.
    /// if --multi is set and -o is a directory, files will be written to outdir as: outdir/chr-start-end.vmatrix.csv.
    /// if --multi is set and -o is a string, the string will be used as a prefix, and files will be written as: <prefix>chr-start-end.vmatrix.csv.
    ///
    /// Examples:
    ///
    /// vplot reads.bam regions.bed > output.vmatrix.csv
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
    ///      - outdir/chr1-1000-2000.vmatrix.csv
    ///     
    ///      - outdir/chr2-1000-2000.vmatrix.csv
    ///
    /// vplot -m -o myPrefix_ reads.bam regions.bed
    ///
    ///   returns:
    ///
    ///      - myPrefix_chr1-1000-2000.vmatrix.csv
    ///
    ///      - myPrefix_chr2-1000-2000.vmatrix.csv
    ///
    /// vplot -m -o outdir/myPrefix_ reads.bam regions.bed
    ///
    ///   returns:
    ///
    ///      - outdir/myPrefix_chr1-1000-2000.vmatrix.csv
    ///
    ///      - outdir/myPrefix_chr2-1000-2000.vmatrix.csv
    ///
    #[structopt(default_value = "-", short = "o", long)]
    output: String,
    /// Write an interactive vplot heatmap in HTML format.
    /// Files are suffixed with `.vplot.html` following the rules outlined in the `--multi` helptext.
    #[structopt(short = "ht", long = "html")]
    write_html: bool,
    //TODO: add flag to use bed column as matrix output name

}

fn connect_indexed_bam(path: &std::path::PathBuf) -> bam::IndexedReader {
    bam::IndexedReader::from_path(path)
        .ok()
        .expect("Cannot read bam file. Ensure an index is present (samtools index <bamfile>).")
}

/// if path is - set prefix to ""
fn get_output_prefix(path: &str) -> &str {
    match path {
        "-" => "",
        _ => path
    }
}

/// If path is a directory, return "{path}/{bam}.vmatrix.csv",
/// otherwise ignore bam name & return only path.
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

/// If path is a directory, return "{path}/{bam}.vplot.html",
/// otherwise ignore bam name & return only "{path}.vplot.html".
fn get_output_html_filename(path: &str, bam: &str) -> String {
    // if path is directory write to dir/vmatrix.csv
    // https://stackoverflow.com/a/64952550

    let outfile;

    if path.ends_with("/") {
        outfile = format!("{}/{}.vplot.html", &path, &bam);
    } else {
        //outfile = path.to_string();
        outfile = format!("{}.vplot.html", &path);
    }

    outfile

}

fn main() {
    let args = Cli::from_args();
    let regions = VRegions::new(&args.regions);
    let bed_reader = bed::Reader::from_file(regions.path);
    let mut bam_reader = connect_indexed_bam(&args.bam);
    let output_prefix = get_output_prefix(&args.output);
    let write_html = args.write_html;

    let mut ventry = VEntry::new();

    let mut vmatrix = VMatrix::new(&regions, &args.max_fragment_size, &args.multi, &args.invert);

    for bed_region in bed_reader.unwrap().records() {
        let region = bed_region.ok()
                               .expect("Error reading bed region");

        //println!("{}\t{}\t{}", region.chrom(), region.start(), region.end());

        let region_chrom = region.chrom();
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
                ventry.update(&region, &read);

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

        // If -multi = True, write matrix to file at "{prefix}{chrom}-{start}-{end}.csv"
        // then zero matrix to hold signal for next region
        if !vmatrix.aggregate {

            // write csv named "{prefix}{chr}-{start}-{end}.csv"
            let region_name = format!("{}{}-{}-{}.vmatrix.csv", output_prefix, region_chrom, region.start(), region.end());

            if write_html {
                let html_path = format!("{}{}-{}-{}.vplot.html", output_prefix, region_chrom, region.start(), region.end());
                //let mut plot = vmatrix.to_heatmap();
                //plot.to_html(&html_name);
                vmatrix.write_heatmap_html(&html_path)

            }

            vmatrix.write_csv(&region_name);
        }

        // Advance to next bed region
    }

    //TODO:
    // test heatmap:
    //let hm = HeatMap::new_z(vmatrix.matrix.genrows().collect().to_vec());
    //let hm = HeatMap::new_z(vmatrix.matrix.outer_iter().collect::<i64>().to_vec());
    //let heat_matrix_vec = vmatrix.matrix.outer_iter().map(|x| x.to_vec()).collect::<Vec<Vec<i64>>>();
    //let hm = HeatMap::new_z(heat_matrix_vec);

    let bampath = &args.bam.into_os_string().into_string().unwrap();

    if write_html {
        let html_path : String = get_output_html_filename(output_prefix, &bampath);
        vmatrix.write_heatmap_html(&html_path);
    }

    let output_path : String = get_output_filename(output_prefix, &bampath);

    if vmatrix.aggregate {
        match output_prefix {
            "" => vmatrix.write_csv_stdout(),
             _ => vmatrix.write_csv(&output_path)
        }
    }
}
