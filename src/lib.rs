extern crate bio;
extern crate csv;
extern crate ndarray;


use std::convert::TryInto;
use ndarray::prelude::*;
use ndarray_csv::Array2Writer;
use rust_htslib::bam;
use bio::io::bed;

use plotly::HeatMap;
use plotly::common::{ColorScale, ColorScalePalette};

pub struct VEntry {
    /// Start position of fragment
    start: i64,
    /// Insert size of the fragment
    insert_size: i64,
    /// Start coordinate of the bed region
    region_start: i64,
    /// End coordinate of the bed region
    region_end: i64,
}

impl VEntry {
    /// Initialize struct with -1's
    pub fn new() -> VEntry {

        VEntry {
            /// Start coordinate of read
            start: -1,
            /// Total fragment length (TLEN)
            insert_size: -1,
            /// Start point of fragment
            region_start: -1,
            /// End point of fragment
            region_end: -1,
        }
    }

    /// Computes the midpoint of the fragment
    pub fn midpoint(&self) -> i64 {
        // NOTE: midpoint is rounded down if insert_size is odd
        // TODO: Fix?? Round up?
        self.start + (self.insert_size/2)
    }

    /// Return column coordinate of read start position
    /// This is really just 0-based relative coordinates from start of the region.
    pub fn start_col(&self) -> i64 {
        self.start - self.region_start
    }

    /// Return column coordinate of read end position
    /// This is really just 0-based relative coordinates from start of the region
    pub fn end_col(&self) -> i64 {
        (self.start + self.insert_size) - self.region_start
    }

    /// Update fragment entry
    /// Allows operation in same allocation
    pub fn update(&mut self, region: &bed::Record, read: &bam::Record) {
        self.start = read.pos();
        self.insert_size = read.insert_size();
        self.region_start = region.start().try_into().unwrap();
        self.region_end = region.end().try_into().unwrap();

    }

    pub fn row(&self) -> usize {
        // Row in the VMatrix this entry should be inserted into
        // NOTE: region count has to be > 0, else will get all 0's for first region row positions
        // But ndarray is 0 indexed, so subtract 1 from final to get row #
        //((self.insert_size.abs()) - 1).try_into().unwrap()
        self.insert_size.abs().try_into().unwrap()
    }
}

/// Holds the vplot matrix
pub struct VMatrix {
    pub matrix: ndarray::Array2<i64>,
    pub ncol: i64,
    //nrow: i64,
    pub max_fragment_size: i64,
    pub aggregate: bool,
    /// Whether entries should be inserted using inverted y-coordinates
    pub invert: bool
}

impl VMatrix {
    pub fn new(regions: &VRegions, max_fragment_size: &i64, multi: &bool, invert: &bool) -> VMatrix {

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

    pub fn insert_midpoint(&mut self, entry: &VEntry) {
        //TODO: consider doing the ncol check as in insert_fragment_ends?
        if entry.insert_size.abs() <= self.max_fragment_size {
            // Column in the VMatrix this entry should be inserted into
            // value must be usize because matrix library needs usize
            let col = entry.midpoint() - entry.region_start;

            let row = if self.invert {entry.row()} else {self.max_fragment_size as usize - entry.row()};

            // Possible valid column values
            let col_range = 0..self.ncol;
            if col_range.contains(&col) {
                //self.matrix[[entry.row(), col as usize]] += 1;
                self.matrix[[row, col as usize]] += 1;
            }

        }
    }

    pub fn insert_fragment_ends(&mut self, entry: &VEntry) {
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

    pub fn insert_fragment(&mut self, entry: &VEntry) {
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

    /// re-zeroes matrix
    pub fn clear(&mut self) {
        self.matrix.iter_mut().for_each(|m| *m = 0)
    }

    // TODO: support ability to export to multiple filetypes.

    /// writes matrix in csv format to stdout
    /// matrix contents stay in memory
    pub fn write_csv_stdout(&self){
        let mut writer = csv::Writer::from_writer(std::io::stdout());
        writer.serialize_array2(&self.matrix)
              .expect("cannot write matrix to stdout");

    }

    /// Writes matrix to csv at path, then zeros the matrix
    pub fn write_csv(&mut self, path: &str){

        let mut writer = csv::Writer::from_path(&path)
            .expect(&format!("cannot open connection to {} for writing", &path));

        writer.serialize_array2(&self.matrix)
              .expect(&format!("cannot write matrix to {}", &path));

        self.clear();

    }

    pub fn to_heatmap(&self) -> plotly::Plot {

        // need to rev because rows are pushed from bottom-up to plot
        let heatmap = HeatMap::new_z(self.matrix
                       .outer_iter()
                       .map(|x| x.to_vec())
                       .rev()
                       .collect::<Vec<Vec<i64>>>())
                       .color_scale(ColorScale::Palette(ColorScalePalette::Greys));

        //heatmap.color_scale(colors);
        // .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
        let mut plot = plotly::Plot::new();
        plot.add_trace(heatmap);
        plot

    }

    pub fn write_heatmap_html(&self, path: &str){
        self.to_heatmap().to_html(&path);
    }

    //if vmatrix.aggregate {
    //    match output_prefix {
    //        "" => write_mat_stdout(&vmatrix),
    //        _  => write_mat_path(&vmatrix, output_prefix, &args.bam.into_os_string().into_string().unwrap())
    //    }
    //}

}


pub struct VRegions<'a> {
    pub path: &'a std::path::PathBuf,
    pub n_regions: i64,
    width: i64,
}

impl VRegions<'_> {
    pub fn new(path: &std::path::PathBuf) -> VRegions {
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
            let width : i64 = (record.end() - record.start()).try_into().unwrap();

            if set_width != width {
                panic!("All bed regions must be equal width\nWidth: {} of bed entry on line {} does not match width of first entry: {}", width, n_regions, set_width);
            }

        }

        (n_regions, set_width)
    }

}







#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
