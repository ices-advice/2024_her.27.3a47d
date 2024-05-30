## Prepare plots for report

## Before: results_sf.RData (model)
## After:  summary.png (report)

rm(list=ls())

library(icesTAF)
taf.library(FLCore)
taf.library(stockassessment)
taf.library(FLSAM)
library(ggplot2)
library(tidyverse)
library(rmarkdown)

mkdir("report")

data.source     <- file.path(".", "bootstrap",'data')
data.save       <- file.path(".", "data")
model.save      <- file.path(".", "model")
output.save     <- file.path(".", "output")

################################################################################
### ============================================================================
### running markdowns
### ============================================================================
################################################################################

prefix <- 'markdown'
mkdir(file.path("report",prefix))

render(input = file.path("report_markdown_20240315_subgroup.Rmd"),output_dir=file.path("report",prefix))
render(input = file.path("report_markdown_20240316_plenary.Rmd"),output_dir=file.path("report",prefix))
render(input = file.path("report_markdown_figures.Rmd"),output_dir=file.path("report",prefix))
render(input = file.path("report_markdown_tables.Rmd"),output_dir=file.path("report",prefix))
