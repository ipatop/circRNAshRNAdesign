## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
#install library
devtools::install_github("ipatop/circRNAshRNAdesign")
library(circRNAshRNAdesign)

## -----------------------------------------------------------------------------
head(read.delim("../test/circs_totest.txt"))

## -----------------------------------------------------------------------------
OligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv")

## -----------------------------------------------------------------------------
oligos<-OligoDesigner(input_coordinates = "../test/circs_totest.txt",output = "../test/New_out.tsv",writetab = F)

## -----------------------------------------------------------------------------
head(oligos)

