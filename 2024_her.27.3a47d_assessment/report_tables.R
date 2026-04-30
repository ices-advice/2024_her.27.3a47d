## Prepare tables for report

## Before: catage.csv, catage_full.csv, survey_heras_full.csv (data),
##         fatage.csv, natage.csv, summary.csv (output)
## After:  catage.csv, catage_full.csv, fatage.csv, natage.csv, summary.csv,
##         survey_heras_full.csv (report)

rm(list=ls())

library(icesTAF)

mkdir("report")

data.save       <- file.path(".", "data")
output.save     <- file.path(".", "output")

### ============================================================================
### Read data
### ============================================================================

## catage (divide, round)
catage <- read.taf(file.path(data.save,"catage.csv"))
catage <- plus(round(div(catage, -1)))

## catage.full (trim year, row sum, divide, round)
catage.full <- read.taf(file.path(data.save,"catage_full.csv"))
catage.full <- plus(tail(catage.full, 16))
catage.full$Total <- rowSums(catage.full)
catage.full <- round(div(catage.full, -1))

## survey.heras.full (row sum, divide, round)
survey.heras.full <- read.taf(file.path(data.save,"survey_heras_full.csv"), na.strings="-1")
survey.heras.full$Total <- rowSums(survey.heras.full[-1])
survey.heras.full <- round(div(survey.heras.full, -1))

### ============================================================================
### Read Outputs
### ============================================================================

## natage (divide, round)
natage <- read.taf(file.path(output.save,"natage_sf.csv"))
natage <- plus(round(div(natage, -1)))

## fatage (round)
fatage <- read.taf(file.path(output.save,"fatage_sf.csv"))
fatage <- plus(rnd(fatage, -1, 3))

## summary (divide, round)
summary <- read.taf(file.path(output.save,"summary_sf.csv"))
summary <- div(summary, "Rec|TSB|SSB|Catch|Landings", grep=TRUE)
summary <- rnd(summary, "Rec|TSB|SSB|Catch|Landings", grep=TRUE)
summary <- rnd(summary, "Fbar|SOP", 3, grep=TRUE)

### ============================================================================
### Export tables
### ============================================================================

setwd("report")
write.taf(catage.full)        # 2.7
write.taf(survey.heras.full)  # 3.1.3
write.taf(catage)   # 6.3.1
write.taf(summary)  # 6.3.12
write.taf(fatage)   # 6.3.13
write.taf(natage)   # 6.3.14
setwd("..")

### ============================================================================
### Read Outputs
### ============================================================================

