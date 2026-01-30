# Markovian-Spatiotemporal-Propagation
This Repository contains the Reproducibility Material of "_Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee). The following includes a roadmap for this repository, which follows the Workflow to reproduce the analyses. Comprehensive descriptions, and suggestions, for performing the analyses are provided subsequently.
In addition, also the functions implemented in the package `spFFBS` ([**Luca Presicce**](https://lucapresicce.github.io/)) are available in the public Gihub repository [spFFBS R package repository](https://github.com/lucapresicce/spFFBS).

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

This section provides an extensive Workflow to reproduce all the numbers, and figures displayed in "_Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling_". The Workflow is presented separately for each Section, and anticipated by a suggested setup to ease the execution of the analyses.

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

### Section 4.1 - Amortized Bayesian Forecast

Running [`ABF-genfun.R`](ABF-genfun.R), [`ABF-simulation.R`](ABF-simulation.R) produce the results, contained in the following objects: 
* _interpolation plots_: [`heatmap_amortized_Om.png`](./plots/heatmap_amortized_Om.png), [`heatmap_amortized_Y.png`](./plots/heatmap_amortized_Y.png).

In this section is displayed [`heatmap_amortized_Om.png`](./plots/heatmap_amortized_Om.png) as a Figure. While we present [`heatmap_amortized_Y.png`](./plots/heatmap_amortized_Y.png) in the Appendix E.

### Section 4.2 - Model class influence

Running [`MclosedMopen-simulation.R`](MclosedMopen-simulation.R), [`MclosedMopen-graphics.R`](MclosedMopen-graphics.R) produce the results, contained in the following objects: 
* _replication results_: `replication_results.RData`;
* _posterior metrics plot_: [`plot_theta.png`](./plot_theta.png), [`plot_omega.png`](./plot_omega.png), [`plot_sigma.png`](./plot_sigma.png);
* _predictive metrics plot_: [`plot_pred.png`](./plots/plot_pred.png).

In this section are displayed [`plot_theta.png`](./plot_theta.png), [`plot_omega.png`](./plot_omega.png) as Figures, and the contents of 100 replications, collected in `replication_results.Rdata`. While we present [`plot_sigma.png`](./plot_sigma.png), and [`plot_pred.png`](./plot_pred.png) in the Appendix E.

**Note:** The output file `replications_results.RData` is **not included in this repository** because its size exceeds GitHub's 100 MB limit (the file is approximately 210 MB). However, it is **fully reproducible** by running the script [`MclosedMopen-simulation.R`](../script/MclosedMopen-simulation.R).  
Please be aware that this script may take a **long time to execute**, depending on your system’s resources. If needed, the original `replications_results.RData` file can be provided upon request.

### Section 5 - Copernicus case study analysis

Running [`Copernicusdata-analysis.R`](Copernicusdata-analysis.R) produces the results, contained in the following objects: 
* _data analysis results_: [`copernicus_temporal_forecast_points.png`](./copernicus_temporal_forecast_points.png);
* _interpolation & uncertainty quantification plots_: [`copernicus_interpolation_outsample_true.png`](./copernicus_interpolation_outsample_true.png), [`copernicus_interpolation_outsample_pred.png`](./copernicus_interpolation_outsample_pred.png), [`copernicus_forecast_temp_true.png`](./copernicus_forecast_temp_true.png),[`copernicus_forecast_temp_pred.png`](./copernicus_forecast_temp_pred.png), [`copernicus_forecast_rain_true.png`](./copernicus_forecast_rain_true.png),[`copernicus_forecast_rain_pred.png`](./copernicus_forecast_rain_pred.png), [`copernicus_forecast_wind_true.png`](./copernicus_forecast_wind_true.png), [`copernicus_forecast_wind_pred.png`](./copernicus_forecast_wind_pred.png), [`copernicus_forecast_evps_true.png`](./copernicus_forecast_evps_true.png), [`copernicus_forecast_evps_pred.png`](./copernicus_forecast_evps_pred.png), [`copernicus_interpolation_insample_true.png`](./copernicus_interpolation_insample_true.png), [`copernicus_interpolation_insample_pred.png`](./copernicus_interpolation_insample_pred.png).

In this section are displayed [`copernicus_temporal_forecast_points.png`](./copernicus_temporal_forecast_points.png), and  [`copernicus_interpolation_outsample_true.png`](./copernicus_interpolation_outsample_true.png), [`copernicus_interpolation_outsample_pred.png`](./copernicus_interpolation_outsample_pred.png) as Figures. While Figures [`copernicus_forecast_temp_true.png`](./copernicus_forecast_temp_true.png),[`copernicus_forecast_temp_pred.png`](./copernicus_forecast_temp_pred.png), [`copernicus_forecast_rain_true.png`](./copernicus_forecast_rain_true.png),[`copernicus_forecast_rain_pred.png`](./copernicus_forecast_rain_pred.png), [`copernicus_forecast_wind_true.png`](./copernicus_forecast_wind_true.png), [`copernicus_forecast_wind_pred.png`](./copernicus_forecast_wind_pred.png), [`copernicus_forecast_evps_true.png`](./copernicus_forecast_evps_true.png), [`copernicus_forecast_evps_pred.png`](./copernicus_forecast_evps_pred.png), [`copernicus_interpolation_insample_true.png`](./copernicus_interpolation_insample_true.png), [`copernicus_interpolation_insample_pred.png`](./copernicus_interpolation_insample_pred.png) are described in the Section 6 body, but presented in the Appendix E.

### Appendix C.1 - Space-time weights dynamics

Running [`Weightsdynamics-simulation.R`](Weightsdynamics-simulation.R) produces the results, contained in the following object: 
* _weights dynamics_: [`plot_weights_dynamic.png`](./plot_weights_dynamic.png);
* _parameter dynamics_: [`plot_par_dynamic.png`](./plot_par_dynamic.png);
* _weights distribution_: [`plot_weight_distr.png`](./plot_weight_distr.png);

In this section are displayed [`plot_weights_dynamic.png`](./plot_weights_dynamic.png), [`plot_weight_distr.png`](./plot_weight_distr.png), [`plot_par_dynamic.png`](./plot_par_dynamic.png) as Figures.

### Appendix D - Copernicus exploratory data analysis

Running [`Copernicusdata-eda.R`](Copernicusdata-eda.R) produces the results, contained in the following objects: 
* _boxplots & histograms_: [`box_resp_temp.png`](./box_resp_temp.png), [`box_resp_rain.png`](./box_resp_rain.png), [`box_resp_wind.png`](./box_resp_wind.png), [`box_resp_evps.png`](./box_resp_evps.png), [`hist_resp_temp.png`](./hist_resp_temp.png), [`hist_resp_rain.png`](./hist_resp_rain.png), [`hist_resp_wind.png`](./hist_resp_wind.png), [`hist_resp_evps.png`](./hist_resp_evps.png);
* _correlation analysis_: [`copernicus_eda_corr.png`](./copernicus_eda_corr.png);
* _spatiotemporal variograms_: [`copernicus_eda_stvariogram.png`](./copernicus_eda_stvariogram.png);
* _Hovmöller diagrams_: [`copernicus_eda_hovmoller_lat.png`](./copernicus_eda_hovmoller_lat.png), [`copernicus_eda_hovmoller_lon.png`](./copernicus_eda_hovmoller_lon.png).

In this section are displayed [`copernicus_eda_corr.png`](./copernicus_eda_corr.png), [`copernicus_eda_stvariogram.png`](./copernicus_eda_stvariogram.png), [`copernicus_eda_hovmoller_lat.png`](./copernicus_eda_hovmoller_lat.png), [`copernicus_eda_hovmoller_lon.png`](./copernicus_eda_hovmoller_lon.png) as Figures, and [`box_resp_temp.png`](./box_resp_temp.png), [`box_resp_rain.png`](./box_resp_rain.png), [`box_resp_wind.png`](./box_resp_wind.png), [`box_resp_evps.png`](./box_resp_evps.png), [`hist_resp_temp.png`](./hist_resp_temp.png), [`hist_resp_rain.png`](./hist_resp_rain.png), [`hist_resp_wind.png`](./hist_resp_wind.png), [`hist_resp_evps.png`](./hist_resp_evps.png) are collected in a Table.

--------------------------------------------------------------------------------
## Contacts

| **Author**|**Maintainer** |**Reference** |
| :--- | :--- | :--- |
| Luca Presicce (l.presicce@campus.unimib.it), Sudipto Banerjee (sudipto@ucla.edu) | Luca Presicce (l.presicce@campus.unimib.it) | "_Adaptive Markovian Spatiotemporal Transfer Learning in Multivariate Bayesian Modeling_" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee)  |



 .
