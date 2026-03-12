# Class 17: Analyzing sequencing data in the cloud
Kyle Canturia (A17502778)

\##Downstream analysis

``` r
#BiocManager::install("tximport")
library(tximport)

# setup the folder and filenames to read
folders <- dir(pattern="SRR21568*")
samples <- sub("_quant", "", folders)
files <- file.path( folders, "abundance.h5" )
names(files) <- samples

#txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
```
