## Conversion script

lfmm2baypass.R is an R script that converts a genotype matrix in the lfmm format, with individuals in rows and genotyped loci in colunms, in a format appropriate for using baypass.

## Installing extendedForest and gradientForest in R

The extendedForest package is required to install gradientForest in R. Thanks to Roman, there is a corrected source for extendedForest that can be installed like this (linux and macos).

Run the following command from the terminal

R CMD INSTALL extendedForest_1.6.2.tar.gz

Then open R or Rstudio and install gradientForest like this

install.packages("gradientForest", repos="<http://R-Forge.R-project.org>")

## Installing gfortran for macOS

it worked with gfortran-12.2-universal.pkg  from <https://cran.r-project.org/bin/macosx/tools/>
