# Markovian-Spatiotemporal-Propagation
This Repository contains the Reproducibility Material of "_Dynamic Bayesian Predictive Stacking via Markovian Spatiotemporal Propagation_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee). The following includes a roadmap for this repository, which follows the Workflow to reproduce the analyses. Comprehensive descriptions and suggestions for performing the analyses are provided subsequently.
In addition, the functions implemented in the package `spFFBS` ([**Luca Presicce**](https://lucapresicce.github.io/)) are available in the public GitHub repository [spFFBS R package repository](https://github.com/lucapresicce/spFFBS).

<!--
Novel approach to performing online learning for multivariate spatiotemporal models. We aim to build a Markovian dependence structure between the incoming datasets at each time instant. In doing so, we exploit the suitable matrix formulation obtainable for dynamic linear models.
-->

--------------------------------------------------------------------------------
## Roadmap of the Repository

| Folder | Contents |
| :--- | :---: |
| **data** | preprocessed dataset in `.Rdata` format & preprocessing scripts |
| **plots** | data analyses/simulations results in `.Rdata` format & figures in paper/supplement  |
| **script** | data analyses/simulations working scripts in `.R` format |


---------------------------------------------------------------------------------------------------------------------------
## Workflow for Reproducible Results

This section provides an extensive Workflow to reproduce all the numbers and figures displayed in "_Dynamic Bayesian Predictive Stacking via Markovian Spatiotemporal Propagation_". The Workflow is presented separately for each Section, and anticipated by a suggested setup to ease the execution of the analyses.

### Working directory

Since the structure of the R Scripts, the computations are organized considering the starting working directory of the entire repository. The scripts begin with:
```{r, echo = F, eval = F, collapse = TRUE}
setwd(".../Markovian-Spatiotemporal-Propagation")
```
where `".../"` represents the path on the user's machine, and then, the directory path where this repository should be placed before executing the scripts. The best option is to clone this repository on the local machine, by executing the following block of code in a `shell`. Once the location to clone this repository is chosen, open the command line and execute:
```{sh}
git clone https://github.com/lucapresicce/Markovian-Spatiotemporal-Propagation.git
```
If not possible, it is possible to execute the scripts by omitting the `setwd("...")` command, but it is mandatory to create the _plots_ folder in the working directory. This allows you to save the results and figures directly inside it.

### Package environments

The most important is the 'spFFBS' package, for which installation of the `devtools` library is required:
```{r}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once devtools is available on the local machine, installation from the GitHub repository proceeds as follows:
```{r}
devtools::install_github("lucapresicce/spFFBS")
```
The [`script`](./script) folder also includes the Rcpp source [`FFBS-DYNBPS-struct-v2.cpp`](./script/FFBS-DYNBPS-struct-v2.cpp) and its R wrapper [`FFBS-DYNBPS-v2.R`](./script/FFBS-DYNBPS-v2.R), compiled automatically at runtime via `Rcpp::sourceCpp()`. These files are provided to ensure full reproducibility ahead of the forthcoming `spFFBS` CRAN new release, which will incorporate these routines natively.

### Section 4.1 - Amortized Bayesian Forecast

Running [`ABF-genfun.R`](./script/ABF-genfun.R), [`ABF-simulation.R`](./script/ABF-simulation.R) produce the results, contained in the following objects: 
* _interpolation plots_: [`heatmap_amortized_Om.png`](./plots/heatmap_amortized_Om.png), [`heatmap_amortized_Y.png`](./plots/heatmap_amortized_Y.png).

In this section, [`heatmap_amortized_Om.png`](./plots/heatmap_amortized_Om.png) is displayed as a Figure. We present [`heatmap_amortized_Y.png`](./plots/heatmap_amortized_Y.png) in Supplementary Section 5.

### Section 5 - Copernicus case study analysis

Running [`Copernicusdata-analysis.R`](./script/Copernicusdata-analysis.R) produces the results, contained in the following objects:
* _data analysis results_: [`copernicus_temporal_forecast_points.png`](./plots/copernicus_temporal_forecast_points.png), [`copernicus_temporal_forecast_lines.png`](./plots/copernicus_temporal_forecast_lines.png);
* _interpolation & uncertainty quantification plots_: [`fig7_IS_spatial.png`](./plots/fig7_IS_spatial.png), [`fig7_OS_spatial.png`](./plots/fig7_OS_spatial.png), [`figS12_temp.png`](./plots/figS12_temp.png), [`figS13_rain.png`](./plots/figS13_rain.png), [`figS14_wind.png`](./plots/figS14_wind.png), [`figS15_evps.png`](./plots/figS15_evps.png).

In this section are displayed [`copernicus_temporal_forecast_points.png`](./plots/copernicus_temporal_forecast_points.png), [`copernicus_temporal_forecast_lines.png`](./plots/copernicus_temporal_forecast_lines.png), and [`fig7_OS_spatial.png`](./plots/fig7_OS_spatial.png) as Figures. While Figures [`fig7_IS_spatial.png`](./plots/fig7_IS_spatial.png), [`figS12_temp.png`](./plots/figS12_temp.png), [`figS13_rain.png`](./plots/figS13_rain.png), [`figS14_wind.png`](./plots/figS14_wind.png), [`figS15_evps.png`](./plots/figS15_evps.png) are presented in the Supplementary Section 5.

**Note:** The input data [`data/copernicus_data.Rdata`](./data/copernicus_data.Rdata) and [`data/copernicus_predictors.Rdata`](./data/copernicus_predictors.Rdata) are available in the repository. The hyperparameter grid [`data/variogram_informed_grid.Rdata`](./data/variogram_informed_grid.Rdata) is produced by the EDA script below but is available directly without being generated first.

### Supplementary Section 3.1 - Model class influence

Running [`MclosedMopen-simulation.R`](./script/MclosedMopen-simulation.R) produces the results and all associated figures, contained in the following objects:
* _replication results_: `replication_results.RData`;
* _posterior metrics plot_: [`plot_theta_C.png`](./plots/plot_theta_C.png), [`plot_theta_O.png`](./plots/plot_theta_O.png), [`plot_omega_C.png`](./plots/plot_omega_C.png), [`plot_omega_O.png`](./plots/plot_omega_O.png), [`plot_sigma.png`](./plots/plot_sigma.png);
* _predictive metrics plot_: [`plot_pred.png`](./plots/plot_pred.png).

In this section are displayed [`plot_theta_C.png`](./plots/plot_theta_C.png), [`plot_theta_O.png`](./plots/plot_theta_O.png), [`plot_omega_C.png`](./plots/plot_omega_C.png), [`plot_omega_O.png`](./plots/plot_omega_O.png), [`plot_sigma.png`](./plots/plot_sigma.png), and [`plot_pred.png`](./plots/plot_pred.png) as Figures, and the contents of 50 replications, collected in `replication_results.Rdata`.

**Note:** The output file `replications_results.RData` is **not included in this repository** because its size exceeds GitHub's 100 MB limit (the file is approximately 9 GB). However, it is **fully reproducible** by running the script [`MclosedMopen-simulation.R`](./script/MclosedMopen-simulation.R). Please be aware that this script may take a **long time to execute**, depending on your system’s resources. If needed, the original `replications_results.RData` file can be provided upon request.

### Supplementary Section 3.2 - Space-time weights dynamics

Running [`Weightsdynamics-simulation.R`](./script/Weightsdynamics-simulation.R) produces the results, contained in the following object:
* _weights dynamics_: [`plot_weights_dynamic.png`](./plots/plot_weights_dynamic.png);
* _parameter dynamics_: [`plot_par_dynamic.png`](./plots/plot_par_dynamic.png);
* _weights distribution_: [`plot_weight_distr.png`](./plots/plot_weight_distr.png);

In this section are displayed [`plot_weights_dynamic.png`](./plots/plot_weights_dynamic.png), [`plot_weight_distr.png`](./plots/plot_weight_distr.png), [`plot_par_dynamic.png`](./plots/plot_par_dynamic.png) as Figures.

### Supplementary Section 4 - Copernicus exploratory data analysis

Running [`Copernicusdata-eda.R`](./script/Copernicusdata-eda.R) produces the results, contained in the following objects:
* _boxplots & histograms_: [`box_resp_01.png`](./plots/box_resp_01.png), [`box_resp_02.png`](./plots/box_resp_02.png), [`box_resp_03.png`](./plots/box_resp_03.png), [`box_resp_04.png`](./plots/box_resp_04.png), [`hist_resp_01.png`](./plots/hist_resp_01.png), [`hist_resp_02.png`](./plots/hist_resp_02.png), [`hist_resp_03.png`](./plots/hist_resp_03.png), [`hist_resp_04.png`](./plots/hist_resp_04.png);
* _correlation analysis_: [`copernicus_eda_corr.png`](./plots/copernicus_eda_corr.png);
* _spatiotemporal variograms_: [`copernicus_eda_stvariogram.png`](./plots/copernicus_eda_stvariogram.png);
* _Hovmöller diagrams_: [`copernicus_eda_hovmoller_lat.png`](./plots/copernicus_eda_hovmoller_lat.png), [`copernicus_eda_hovmoller_lon.png`](./plots/copernicus_eda_hovmoller_lon.png).

In this section are displayed [`copernicus_eda_corr.png`](./plots/copernicus_eda_corr.png), [`copernicus_eda_stvariogram.png`](./plots/copernicus_eda_stvariogram.png), [`copernicus_eda_hovmoller_lat.png`](./plots/copernicus_eda_hovmoller_lat.png), [`copernicus_eda_hovmoller_lon.png`](./plots/copernicus_eda_hovmoller_lon.png) as Figures, and [`box_resp_01.png`](./plots/box_resp_01.png)–[`box_resp_04.png`](./plots/box_resp_04.png), [`hist_resp_01.png`](./plots/hist_resp_01.png)–[`hist_resp_04.png`](./plots/hist_resp_04.png) are collected in a Table (variables ordered as: Temperature, Rain, Wind, Evapotranspiration).

This script also produces [`data/variogram_informed_grid.Rdata`](./data/variogram_informed_grid.Rdata) and should be run **before** the data analysis script in Section 5. The input data [`data/copernicus_data.Rdata`](./data/copernicus_data.Rdata) and [`data/copernicus_predictors.Rdata`](./data/copernicus_predictors.Rdata) are available in the repository.

--------------------------------------------------------------------------------
## Contacts

| **Author**|**Maintainer** |**Reference** |
| :--- | :--- | :--- |
| Luca Presicce (l.presicce@campus.unimib.it), Sudipto Banerjee (sudipto@ucla.edu) | Luca Presicce (l.presicce@campus.unimib.it) | "_Dynamic Bayesian Predictive Stacking via Markovian Spatiotemporal Propagation_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee)  |



 .
