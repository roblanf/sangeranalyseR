<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Travis build status](https://travis-ci.org/roblanf/sangeranalyseR.svg?branch=master)](https://travis-ci.org/roblanf/sangeranalyseR) ![documentation build status](https://readthedocs.org/projects/pip/badge/)
<!-- badges: end -->
  
# sangeranalyseR

sangeranalseR is an R package that provides fast, flexible, and reproducible workflows for assembling your sanger seuqencing data into contigs. It adds to a list of already widely-used tools, like Geneious, CodonCode Aligner and Phred-Phrap-Consed. What makes it different from these tools is that it’s free, it’s open source, and it’s in R.

For more information, please check our <b>📒<a href="https://sangeranalyser.readthedocs.io/en/latest/">sangeranalyseR Documentation</a></b>.

<br>

## Citation
sangeranalyseR is on [***Genome Biology and Evolution (GBE)***](https://academic.oup.com/gbe/advance-article/doi/10.1093/gbe/evab028/6137837?guestAccessKey=a28b32d6-ffab-41f2-8132-9c2dd28b99fe) and [***Bioconductor 3.13***](https://bioconductor.org/packages/release/bioc/html/sangeranalyseR.html) now.

If you use sangeranalyseR in your published work, please cite

> **Kuan-Hao Chao, Kirston Barton, Sarah Palmer, and Robert Lanfear (2021). "sangeranalyseR: simple and interactive processing of Sanger sequencing data in R" in Genome Biology and Evolution. DOI: [doi.org/10.1093/gbe/evab028](https://doi.org/10.1093/gbe/evab028)**

<br>

## Quick Start Guide
### 1. Installation

#### (1) System requirements
* R >= 4.0.0 (current)
* Rstudio (recommended)

#### (2) Install from Bioconductor
To install this package, start R (version “4.0”) and enter:
``` R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("sangeranalyseR")
```

#### (3) Install from GitHub
If you haven’t installed the `devtools` package before, please install it first:

``` R
install.packages("devtools")
```

Then run the following code in your R console to install the newest version from Github.

``` R
library(devtools)

# Download it from the master branch
install_github("roblanf/sangeranalyseR", ref = "master")

# Download it from the develop branch
install_github("roblanf/sangeranalyseR", ref = "develop")

library(sangeranalyseR)
```

<br>

### 2. A Reproducible Example

Here we demonstrate a simple and reproducible example for using sangeranalyseR to generate a consensus read from 8 sanger ab1 files (4 contigs and each includes a forward and a reverse read).

#### (1) Prepare your input files & loading

The data of this example is in the sangeranalyseR package; thus, you can simply get its path from the library.

``` R
rawDataDir <- system.file("extdata", package = "sangeranalyseR")
parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica', 'ACHLO')
```

#### (2) Load and analyse your data

Run the following on-liner to create the *SangerAlignment* object.

``` R
ACHLO_contigs <- SangerAlignment(ABIF_Directory     = parentDir,
                                 REGEX_SuffixForward = "_[0-9]*_F.ab1$",
                                 REGEX_SuffixReverse = "_[0-9]*_R.ab1$")
```

#### (3) Explore your data

Launch the Shiny app to check the visualized results.

``` R
launchApp(ACHLO_contigs)
```
The following figure shows the SangerAlignment Shiny app user interface.![image](https://user-images.githubusercontent.com/26069376/126213907-44180b80-29b2-4b6c-9cd3-41f01733cc97.png)

<img src="https://i.imgur.com/gwY6AqB.png" style="width:100%">

#### (4) Output your aligned contigs

Write each contig and the aligned consensus read into FASTA files.

``` R
writeFasta(ACHLO_contigs)
```

#### (5) Generate an interactive report

Last but not least, generate an Rmarkdown report to store all the sequence information.

``` R
generateReport(ACHLO_contigs)
```

