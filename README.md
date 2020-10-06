# A commandline tool for making vplot matrices

## Install

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
Note that in the returned matrix, the columns are genomic coordinate (leftmost column is 5' end), and rows are fragment size (top row is *smallest* fragment size)
``` sh
vplot <bam file> <bed file> --max-size <max fragment size to use>
```
