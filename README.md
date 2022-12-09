
# circRNAshRNAdesign

<!-- badges: start -->
<!-- badges: end -->

The goal of circRNAshRNAdesign is to create the oligos to knockdown specific circRNAs  

Contents
========

 * [Why?](#why)
 * [Installation](#installation)
 * [Usage](#usage)
 
### Why?

I wanted a tool that allows you to:

+ Get the specific _circRNA junstion_.
+ Get the putative _circRNA sequence_.
+ Design shRNAs against circRNAs.
+ Design shRNAs shifts.


### Installation
---

You can install the development version of circRNAshRNAdesign from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ipatop/circRNAshRNAdesign")
```

### Usage

This is a basic example which shows you how to use it

``` r
library(circRNAshRNAdesign)
## basic example code
OligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv")
```
Input should look like this
```r
head(read.delim("../test/circs_totest.txt"))
```

Run to create an output table
```r
OligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv")
```

Run to create a DataFrame, if writetab = F
```r
oligos<-OligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv",writetab = F)
```

Output
```r
head(oligos)
```

