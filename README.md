# Data-based clustering in prediction of cervical cancer DNA methylation using pan-cancer genetic and clinical data

This repository contains code for the paper "Data-based clustering in prediction of cervical cancer DNA methylation using pan-cancer genetic and clinical data" (Pai & Rao, 2025+).

## Constants, libraries, and functions

1. `preliminaries.R`: Load libraries, create constants
2. `cmmp_functions.R`: Contains the code for CMMP (Jiang et. al., 2018) for versions of R that don't have the package `cmmp`, a wrapper function
3. `sim_functions.R`: Generate data, 3 methods (regression prediction, k-means + CMMP, naive estimate), error metrics

## Code for the TCGA data analysis

1. `download_data.R`: Download data via TCGA retriever (takes time)
2. `merge_data.R`: Merge cancer and data types
3. `preprocessing.R`: Filter variables, scale variables, exploratory plots
4. `clustering.R`: Find number of clusters via gap stat, cluster data, plots
5. `modeling.R`: Run regression prediction and k-means + CMMP on the TCGA data
6. `analysis_real.R`: Compile TCGA data results, results plots

## Code for the simulation study

1. `sim_config.R`: Create a list of simulation settings and seeds for a shell script to refer to
2. `sim.R`: Run simulation specified rows of sim_config via a shell script
3. `analysis.R`: Compile results from individual simulations and generate plots
