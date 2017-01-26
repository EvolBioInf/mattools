# mattools — utilities for distance matrices

Phylip distance matrices are a common intermediate format in phylogeny reconstruction. This tool set consists of a small set of utilities for their manipulation, formatting and comparison.

This readme covers all necessary instructions to get the mattools up and running.

## Usage

Let us assume you have `mat` installed via the way described below. Then you get help to `mat` and all its subcommands via the option `--help`.

    $ mat --help
    The available commands are:
     compare     Compute the distance between two matrices
     format      Format the distance matrix
     grep        Print submatrix for names matching a pattern
     nj          Convert to a tree by neighbor joining

    Use 'mat <command> --help' to get guidance on the usage of a command.

### Comparison

To check whether two distance matrices are equal, use `mat compare`. The input matrices will be interpreted as two vectors and their euclidean distance computed. Thus a distance of zero indicates equality. To circumvent problem with differently sized matrices, only those values are included in the computation, whose corresponding names equal.

### Formatting

Unfortunately, the phylip distance matrix format is poorly designed, described, and badly implemented in different tools. With `mat format` all these formatting differences can be removed.

    $ cat lowertriangle.mat
    2
    A  
    B  0.5
    $ mat format lowertriangle.mat
    2
    A          0.0000e+00 5.0000e-01
    B          5.0000e-01 0.0000e+00

To verify that the distance matrix is indeed a distance matrix in the mathematical sense, the option `--validate` can be used. The mat tools will then hunt for errors and try to fix them, where possible.

### Filtering

To remove or extract individual lines and submatrices `mat grep` can be used. It takes a regular expression and checks it against the names. Names, not matching the pattern are discarded from the output. This behaviour can be changed with the flag `--invert-match`.

### Neighbor Join

The mattools also come with a module for building a phylogeny via neighbor joining. The resulting tree also contains support values computed via quartet analysis. See the following paper for a description of the process: [Klötzl & Haubold (2016)](http://www.mdpi.com/2075-1729/6/1/11/htm).


## Building

Download the source and then execute the following steps.

    $ autoreconf -fi
    $ ./configure
    $ make
    $ make install  # try sudo

## License

Copyright © 2017 Fabian Klötzl  
License GPLv3+: GNU GPL version 3 or later.

This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. The full license text is available at <http://gnu.org/licenses/gpl.html>.

Individual files may be licensed differently.

## Contact

In case of bugs or unexpected errors don't hesitate to send me a mail: kloetzl@evolbio.mpg.de
