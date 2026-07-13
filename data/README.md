# Data Availability for Reproducibility 

"_Dynamic Bayesian Predictive Stacking via Markovian Spatiotemporal Propagation_" considers a case study strictly related to critical problems in today's Climate Sciences. As a reminder, the analyses work on multiple atmospheric components, providing a full Bayesian inference over Europe's surface, including millions of units. Raw data is freely available for download from the following web portals: [Copernicus Climate data Store](https://cds.climate.copernicus.eu/), for which metadata is available.

## Preprocessing

Since the original format and size of raw datasets, some preprocessing steps are mandatory before starting any analysis. [`Copernicusdata-cleaning.R`](Copernicusdata-cleaning.R) is the R script (`.R` format) that concerns the preprocessing of Climate data, passing from raw data in `.nc` format (downloadable at the provided link) to the `.Rdata` object.

However, due to the massive dimensions, raw datasets are not loaded in this folder, but the preprocessed data are then available. In order to use this preprocessing script, first raw data must be downloaded from the link.

## Case Study Datasets

Introducing now [`copernicus_data.Rdata`](copernicus_data.Rdata) and [`copernicus_predictors.Rdata`](copernicus_predictors.Rdata), the `.Rdata` objects containing the preprocessed set of outcomes and predictors, respectively, used for the analyses.[`variogram_informed_grid.Rdata`](variogram_informed_grid.Rdata) contains the spatial hyperparameter grid derived from spatiotemporal variogram analysis (see Supplementary Section 4). It is produced by [`Copernicusdata-eda.R`](../script/Copernicusdata-eda.R) and loaded by [`Copernicusdata-analysis.R`](../script/Copernicusdata-analysis.R). It is available directly in this folder without the need to regenerate it. See both Section 5 and Supplementary Section 4 for further details.

The files here are analysis-ready; it suffices to load them (following the related scripts) to perform the analyses, as explained in the Workflow on the main `README.md`.

