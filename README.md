
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scITD Analysis

This repository is used for reproducing all analysis included in our
paper describing the scITD package. See the scITD package here for more
details: <https://github.com/kharchenkolab/scITD>

We use the renv package to ensure that the exact same packages are used
each time the results are reproduced. To get started, you will first
need to install renv.

``` r
devtools::install_version("renv", version = "0.14.0", repos = "http://cran.us.r-project.org")
```

Next, clone this repository and open R in the R project in this folder.
This can be done by clicking on the scITD-Analysis.Rproj file in this
directory in RStudio.

Finally, you can install the packages used in this analysis using the following. This
will store the updated packages in the project library, not affecting your user library.

``` r
renv::restore()
```

When data is loaded in each analysis script, we indicate from which preprocessing
script the processed data was derived from. You will need to download the data from the
source indicated in the respective preprocessing script sources and run the script. 
Every plot from the paper is labeled in the script where it is produced.