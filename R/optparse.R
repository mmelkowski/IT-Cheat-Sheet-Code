#### Parsing command line example from documentation ####

library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true", 
                default=TRUE, help="Print extra output [default]")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false", 
                    dest="verbose", help="Print little output")
parser <- add_option(parser, c("-c", "--count"), type="integer", default=5, 
                help="Number of random normals to generate [default %default]",
                metavar="number")
parse_args(parser, args = c("--quietly", "--count=15"))

## $help
## [1] FALSE
## 
## $verbose
## [1] FALSE
## 
## $count
## [1] 15

# Note that the args argument of parse_args default is commandArgs(trailing=TRUE) so it typically doesn’t need to be explicitly set if writing an Rscript.
# One can also equivalently make options in a list:

library("optparse")
option_list <- list( 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false", 
        dest="verbose", help="Print little output"),
    make_option(c("-c", "--count"), type="integer", default=5, 
        help="Number of random normals to generate [default %default]",
        metavar="number")
    )

parse_args(OptionParser(option_list=option_list), args = c("--verbose", "--count=11"))

## $verbose
## [1] TRUE
## 
## $count
## [1] 11
## 
## $help
## [1] FALSE