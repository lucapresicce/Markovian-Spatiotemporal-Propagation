# Execution Scripts for Reproducibility 

The R scripts collected in this folder allow the entire reproduction of all the numbers and figures shown in "_Dynamic Bayesian Predictive Stacking via Markovian Spatiotemporal Propagation_". Here, a brief reminder of which script is needed to reproduce the contents of a specific Section is presented to ease the reproducibility. Reminding the Workflow on the main `README.md` of this repository, for an exhaustive and detailed approach to reproduce the results.

* Section 4.1: [`ABF-genfun.R`](ABF-genfun.R), [`ABF-simulation.R`](ABF-simulation.R);
* Section 5: [`Copernicusdata-analysis.R`](Copernicusdata-analysis.R);
* Supplementary Section 3.1: [`MclosedMopen-simulation.R`](MclosedMopen-simulation.R);
* Supplementary Section 3.2: [`Weightsdynamics-simulation.R`](Weightsdynamics-simulation.R);
* Supplementary Section 4: [`Copernicusdata-eda.R`](Copernicusdata-eda.R);
* Supplementary Section 5: figures collected from [`ABF-simulation.R`](ABF-simulation.R) and [`Copernicusdata-analysis.R`](Copernicusdata-analysis.R), already listed above.

**Note:** The folder also includes the Rcpp source [`FFBS-DYNBPS-struct-v2.cpp`](FFBS-DYNBPS-struct-v2.cpp) and its R wrapper [`FFBS-DYNBPS-v2.R`](FFBS-DYNBPS-v2.R), compiled automatically at runtime via `Rcpp::sourceCpp()`. These files are provided to ensure full reproducibility ahead of the forthcoming `spFFBS` CRAN new release, which will incorporate these routines natively.
