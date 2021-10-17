# A commandline tool for making V-plot matrices from genomics data

`vplot` is a small Rust tool for generating V-plot matrices. By default it makes
matrices using the original style of fragment midpoints developed in [Henikoff,
et. al. 2011](https://www.pnas.org/content/108/45/18318), but it also supports
making plots from fragment ends, or mapping the whole fragment instead. These
features make `vplot` a versatile tool for visualizing MNase-seq, ATAC-seq, or
CUT&RUN data.

## Install

### Prebuilt binary

Download a compatible binary for your operating system from the [Releases
Page](https://github.com/snystrom/vplot/releases). Can't find one that works?
File an [issue](https://github.com/snystrom/vplot/issues), or compile from
source (instructions below).

Place binary in convenient location (ie ~/bin), add that location to your PATH
if it isn't already.

``` sh
# Adds ~/bin to PATH
echo "export PATH=\$PATH:~/bin" >> ~/.bashrc
```

Test install:

``` sh
vplot --version
```

### Compile from source

[Install Rust](https://www.rust-lang.org/tools/install). It's super easy!

Install the binary
``` sh
cargo install --git https://github.com/snystrom/vplot
```

Test install
``` sh
vplot --version
```

## Usage

Note that in the returned matrix, the columns are genomic coordinate (leftmost
column is 5' end), and rows are fragment size (top row is *smallest* fragment
size)

``` sh
vplot 0.3.11

USAGE:
    vplot [FLAGS] [OPTIONS] <bam> <regions>

FLAGS:
    -h, --help       Prints help information
    -m, --multi      Instead of aggregating reads into 1 matrix, write 1 matrix for each region. Matrices are written as
                     1 csv per region named: `chr-start-end.csv`
    -V, --version    Prints version information

OPTIONS:
    -f, --fragment-type <fragment-type>    How reads are counted in the matrix. Using either the midpoint of the
                                           fragment, fragment ends, or the whole fragment [default: midpoint]  [possible
                                           values: midpoint, ends, fragment]
    -x, --max-size <max-fragment-size>     Maximum fragment size to include in the V-plot matrix [default: 700]
    -o, --output <output>                  Set output file name or output directory. This option behaves differently
                                           depending on which input flags are set. See --help for details [default: -]

ARGS:
    <bam>        Path to an indexed bam file
    <regions>    Path to a bed file (must be in bed4 format: chr, start, end, strand) Of a region (or regions) in
                 which to generate the vplot. If using multiple regions, all entries must be the same width. (NOTE:
                 using more than 1 region is poorly supported currently, but technically works)
```
