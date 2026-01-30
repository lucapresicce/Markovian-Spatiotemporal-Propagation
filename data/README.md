# Data Availability for Reproducibility 

"_Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling_" considers a case study strictly related to critical problems in today's Climate Sciences. As a reminder, the analyses work on multiple atmospheric components, providing a full Bayesian inference over Europe's surface, including millions of units. Raw data is freely available for download from the following web portals: [Copernicus Climate data Store](https://cds.climate.copernicus.eu/), for which metadata is available.

## Preprocessing

Since the original format and size of raw datasets, some preprocessing steps are mandatory before starting any analysis. [`Copernicusdata-cleaning.R`](./data/Copernicusdata-cleaning.R) is the R script (`.R` format) that concerns the preprocessing of Climate data, passing from raw data in `.nc` format (downloadable at the provided link) to the `.Rdata` object;

However, due to the massive dimensions, raw datasets are not loaded in this folder, but the preprocessed data are then available. In order to use this preprocessing script, first raw data must be downloaded from the link.

## Case Study Datasets

Introducing now [`copernicus_data.RData`](./data/copernicus_data.RData), and [`copernicus_predictors.RData`](./data/copernicus_predictors.RData), the `.Rdata` objects containing the preprocessed set of outcomes and predictors, respectively, used for the analyses. See both Section 5 and Appendix D for further details. 

The files here are analysis-ready; it suffices to load them (following the related scripts) to perform the analyses, as explained in the Workflow on the main `README.md`.
