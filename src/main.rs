extern crate bio;

use std::convert::TryInto;
use rust_htslib::{bam, bam::Read};
//use rust_htslib::bam;
use structopt::StructOpt;
use bio::io::bed;

#[derive(StructOpt, Debug)]
#[structopt(setting = structopt::clap::AppSettings::ArgRequiredElseHelp)]
struct Cli {
    #[structopt(parse(from_os_str))]
    bam: std::path::PathBuf,
    #[structopt(parse(from_os_str))]
    regions: std::path::PathBuf,

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
    fn row(&self) -> i64 {
        // NOTE: region count has to be > 0, else will get all 0's for first region row positions
        (self.region_n + 1) * self.insert_size.abs()
    }

    /// Column in the VMatrix this entry should be inserted into
    fn col(&self) -> i64 {
        self.region_end - self.midpoint()
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

fn main() {
    let args = Cli::from_args();
    let bed_reader = bed::Reader::from_file(&args.regions);
    let mut bam_reader = bam::IndexedReader::from_path(&args.bam)
        .ok()
        .expect("Cannot read bam file");

    let mut ventry = VEntry::allocate();

    let mut bed_region_n = 0 as i64;
    for bed_region in bed_reader.unwrap().records() {
        let region = bed_region.ok()
                               .expect("Error reading bed region");

        //println!("{}\t{}\t{}", region.chrom(), region.start(), region.end());

        let region_chrom = region.chrom();
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
                // TODO: write to matrix
                //
                ventry.update(&region, bed_region_n, &record);
                println!("{}\t{}\t{}", ventry.row(), ventry.col(), ventry.midpoint());
            }

            // Advance to next bam region
            check = bam_reader.read(&mut record);
        }

        bed_region_n += 1;
    }

}
