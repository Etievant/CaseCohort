Replication of the simulation studies in "Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data", by Etiévant and Gail (Lifetime Data Analysis, 2024). The methods described in the article are implemented in the file `help.functions.R`.

### Required packages 

```
dplyr, ggplot2, grid, gridExtra, gtable, parallel, survival, xtable.
```

### Scripts

* Script `simulations.R` replicates the simulations proposed by Etiévant and Gail (Lifetime Data Analysis, 2024) in Section 7 and Web Appendix D.2.

* Scripts `simulations_postratifiedweights.R` replicates the simulations in Web Appendix D.5.

* Scripts `simulations_weakerproxy.R` replicates the simulations in Web Appendix D.6.

* Scripts `simulations_missingdata.R` replicates the simulations Web Appendix G.

Each script relies on functions provided in `help.functions.R`.


### Instructions to run each script

* Save the chosen script(s) and file `help.functions.R` in the same directory.

* Open and run the whole script(s).

* The results of the simulations are saved in figures and tables. For example, when running script `simulations.R`, file `details.beta1.csv` will contain the simulation results displayed in Web Tables 3, 9 and 10 in Web Appendix D.2, and file `SimulationResults.RH.Scenario3.csv` will contain the results displayed in Table 1 in Section 7.


### Functions provided in `help.functions.R`

* **estimation** - estimation of log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\left[ \tau_1,\tau_2 \right]$ and pure risk in $\left[ \tau_1,\tau_2 \right]$ and for a given covariate profile $\boldsymbol{x}$.

* **estimation.CumBH** - estimation of log-relative hazard, baseline hazards at each unique event time and cumulative baseline hazard in a given time interval $\left[ \tau_1,\tau_2 \right]$.

* **estimation.PR** - estimation of the pure risk in a given time interval $\left[ \tau_1,\tau_2 \right]$ and for a given covariate profile $\boldsymbol{x}$, from the values of the log-relative hazard and cumulative baseline hazard in $\left[ \tau_1,\tau_2 \right]$.

* **influences** - estimation of influences on log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\left[ \tau_1,\tau_2 \right]$ and on the pure risk in $\left[ \tau_1,\tau_2 \right]$ and for a given covariate profile $\boldsymbol{x}$. Also provides parameters estimates. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.RH** - estimation of influences on log-relative hazard. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.CumBH** - estimation of influences on log-relative hazard, baseline hazards at each unique event time and cumulative baseline hazard in a given time interval $\left[ \tau_1,\tau_2 \right]$. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.PR** - estimation of influences on the on pure risk in a given time interval $\left[ \tau_1,\tau_2 \right]$ and for a given covariate profile $\boldsymbol{x}$, from that on the log-relative hazard and cumulative baseline hazard in $\left[ \tau_1,\tau_2 \right]$. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.missingdata** - estimation of influences on log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\left[ \tau_1,\tau_2 \right]$ and the on pure risk in $\left[ \tau_1,\tau_2 \right]$ and for a given covariate profile $\boldsymbol{x}$, when covariate data is missing for individuals in the phase-two sample. Also provides parameters estimates.

* **influences.RH.missingdata** - estimation of influences on log-relative hazard, when covariate data is missing for individuals in the phase-two sample.

* **influences.CumBH.missingdata** - estimation of influences on log-relative hazard, baseline hazards at each unique event time, and cumulative baseline hazard in a given time interval $\left[ \tau_1,\tau_2 \right]$, when covariate data is missing for individuals in the phase-two sample.

* **influences.PR.missingdata** - estimation of influences on the on pure risk in a given time interval $\left[ \tau_1,\tau_2 \right]$ and for a given covariate profile $\boldsymbol{x}$, from that of the log-relative hazard and cumulative baseline hazard $\left[ \tau_1,\tau_2 \right]$, when covariate data is missing for individuals in the phase-two sample. 

* **robustvariance** - estimation of the robust variance estimate, for parameters such as log-relative hazard, cumulative baseline hazard or covariate specific pure risk. Works with design weights or calibrated weights, or when covariate data is missing for individuals in the phase-two sample.

* **variance** - estimation of the influence-based variance, that follows the complete variance decomposition, for parameters such as log-relative hazard, cumulative baseline hazard or covariate specific pure risk. Works with design weights or calibrated weights.

* **variance.missingdata** - estimation of the influence-based variance, that follows the complete variance decomposition, for parameters such as log-relative hazard, cumulative baseline hazard or covariate specific pure risk, when covariate data is missing for individuals in the phase-two sample.

* **auxiliary.construction** - construction of the auxiliary variables proposed by Breslow et al. (Stat. Biosci., 2009) and Breslow et al. (IMS, 2013). 

* **calibration** - calibration of the design weights using the raking procedure. The Newton Raphson method is used for the optimization.

* **product.covar.weight** - construction of the product of joint design weights and joint sampling indicators covariances, needed for the phase-two component of the variance (with design or calibrated weights).

* **cohort.generation** - generation of a cohort following the design used in the simulation study in Section 7 and Web Appendix E.

* **conf.interval** - construction of 95% confidence interval with parameter estimate and its variance estimate, assuming normality.

* **densityX1X3** and **densityX1** - used for the computation of the constant baseline hazard in the simulations in Section 7 and Web Appendix E.

* **E.RH** and **E.RH.w** - computation of $\text{E} \lbrace \text{exp}\(\beta_1 X_1 + \beta_2 X_2 + \beta_3 X_3 \)\rbrace$ and $\text{E} \lbrace\text{exp}\(\beta_1 X_1 + \beta_2 X_2 + \beta_3 X_3 | \boldsymbol{W} = \boldsymbol{w} \)\rbrace$, for $\boldsymbol{w}$ = {0, 1, 2, 3}, needed to compute the constant baseline hazard and the sampling fractions to be used for an un-stratified and stratified sampling of the case-cohort in the simulations in Section 7 and Web Appendix E.

* **scenarios** - creation of a dataframe with the different scenarios to be investigated in the simulations in Section 7 and Web Appendix E.

* **scenarios.missingdata** - creation of a dataframe with the different scenarios to be investigated in the simulations in Web Appendix G (when covariate data is missing for individuals in the phase-two sample).
  
* **X.generation** - generation of a gaussian variable, with unit variance and given mean.
 
* **detailed.results** - analysis of the simulation results and creation of .csv and .Rdata files to save the details as in Section 7 and Web Appendix E.


### References

Etiévant L, Gail MH (2024) Cox model inference for relative hazard and pure risk from stratified weight-calibrated case-cohort data. Lifetime Data Anal 30, 572–599. https://doi.org/10.1007/s10985-024-09621-2

Breslow NE, Lumley T, Ballantyne CM, Chambless LE, Kulich M (2009a) Improved Horvitz-Thompson
estimation of model parameters from two-phase stratifed samples: applications in epidemiology.
Stat Biosci 1(1):32–49. https://doi.org/10.1007/s12561-009-9001-6

Breslow NE, Lumley T (2013) Semiparametric models and two-phase samples: applications to Cox
regression. From probability to statistics and back: high-dimensional models and processes—a Festschrift in honor of Jon A. Wellner 9:65–78. https://doi.org/10.1214/12-IMSCOLL906


