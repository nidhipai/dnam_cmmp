library(cluster)
library(MASS)
library(tidyverse)
library(ggplot2)
#library(cmmp) # Doesn't exist in MSI R version
library(nlme) # Needed for CMMP
library(viridis)
library(patchwork)

source("code/cmmp_functions.R")

# Map race abbreviations to more common ones
race_abbrev <- c(
  AIAN = "AI/AN", Asian = "Asian", BAA = "Black", White = "White", NHPI = "NHPI"
)

cesc_rates <- c(AIAN = 0.014, Asian = .055, BAA = .246, White = .681, NHPI = .002)

select <- dplyr::select
