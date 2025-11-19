# Data Availability for Reproducibility 

"_Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling_" considers a case study strictly related to critical problems in nowadays Climate Sciences. As a reminder, the analyses work on multiple atmospheric components, providing a full Bayesian inference over Europe's surface, including millions of units. Raw data is freely available for download from the following web portals: [Copernicus Climate data Store](https://cds.climate.copernicus.eu/), for whose metadata is available.

## Preprocessing

Since the original format and size of raw datasets, some preprocessing steps are mandatory before starting any analysis. Here the R script (`.R` format) [`Copernicus data - cleaning.R`](../data/Copernicus data - cleaning.R) concerns the preprocessing of Climate data, passing from raw data in `.nc` format (downloadable at the provided link) to the `.Rdata` object;

However, due to the massive dimensions, raw datasets are not loaded in this folder, but the preprocessed data are then available. In order to use this preprocessing script, first raw data must be downloaded from the link.

## Case Study Datasets

Introducing now [`copernicus_data.RData`](../data/copernicus_data_2024_05_12.RData), the `.Rdata` object containing the preprocessed set of data used for the analyses. See Section 2 for further details. 

The files here are analysis-ready; it suffices to load them (following the related scripts) to perform the analyses, as explained in the Workflow on the main `README.md`.
