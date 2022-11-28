

The attached scripts replicate the simulation studies in "Inference for relative hazard and covariate-specific pure risk estimated from stratified and unstratified case-cohort data", by Etiévant and Gail (2022). The methods described in the article are implemented in the file `help.functions.R`.

### Required packages 

```
dplyr, ggplot2, grid, gridExtra, gtable, parallel, survival, xtable.
```

### Scripts

* Script `simulations.R` replicates the simulations proposed by Etiévant Gail in Section 5 and Web Appendix E.2.

* Scripts `simulations_postratifiedweights.R` replicates the simulations in Web Appendix E.4.

* Scripts `simulations_weakerproxy.R` replicates the simulations in Web Appendix E.5.

* Scripts `simulations_missingdata.R` replicates the simulations Web Appendix I.

Each script relies on functions provided in `help.functions.R`.


### Instructions to run each script

* Save the chosen script(s) and files `help.functions.R` in the same directory.

* Open and run the whole script(s).

* The results of the simulations are saved in figures and tables. For example, when running script `simulations.R`, file details.beta1.csv will contain the simulation results displayed in Web Tables 3, 9 and 15 in Web Appendix E.2, and file SimulationResults.RH.Scenario3.csv will contain the results displayed in Table 1 in Section 5.


### Functions provided in `help.functions.R`

* **estimation** - estimation of log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\[\tau_1,\tau_2\]$ and pure risk in $\[\tau_1,\tau_2\]$ and for a given covariate profile $\boldsymbol{x}$.

* **estimation.CumBH** - estimation of log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\[\tau_1,\tau_2\]$.

* **estimation.PR** - estimation of pure risk in a given time interval $\[\tau_1,\tau_2\]$ and for a given covariate profile $\boldsymbol{x}$, from the values of the log-relative hazard and cumulative baseline hazard in $\[\tau_1,\tau_2\]$.

* **influences** - estimation of influences on log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\[\tau_1,\tau_2\]$ and the on pure risk in $\[\tau_1,\tau_2\]$ and for a given covariate profile $\boldsymbol{x}$. Also provides parameters estimates. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.RH** - estimation of influences on log-relative hazard. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.CumBH** - estimation of influences on log-relative hazard, baseline hazards at each unique event time and cumulative baseline hazard in a given time interval $\[\tau_1,\tau_2\]$. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.PR** - estimation of influences on the on pure risk in a given time interval $\[\tau_1,\tau_2\]$ and for a given covariate profile $\boldsymbol{x}$, from that of the log-relative hazard and cumulative baseline hazard in $\[\tau_1,\tau_2\]$. Can take calibration of the design weight into account if auxiliary variables observations are provided.

* **influences.missingdata** - estimation of influences on log-relative hazard, baseline hazards at each unique event time, cumulative baseline hazard in a given time interval $\[\tau_1,\tau_2\]$ and the on pure risk in $\[\tau_1,\tau_2\]$ and for a given covariate profile $\boldsymbol{x}$, when covariate data is missing for individuals in the case cohort. Also provides parameters estimates.

* **influences.RH.missingdata** - estimation of influences on log-relative hazard.

* **influences.CumBH.missingdata** - estimation of influences on log-relative hazard, baseline hazards at each unique event time, and cumulative baseline hazard in a given time interval $\[\tau_1,\tau_2\]$.

* **influences.PR.missingdata** - estimation of influences on the on pure risk in a given time interval $\[\tau_1,\tau_2\]$ and for a given covariate profile $\boldsymbol{x}$, from that of the log-relative hazard and cumulative baseline hazard $\[\tau_1,\tau_2\]$, when covariate data is missing for individuals in the case cohort. 

* **robustvariance** - estimation of the robust variance estimate, for parameters such as log-relative hazard, cumulative baseline hazard or covariate specific pure-risk. Works with design weights or calibrated weights, or when covariate data is missing for individuals in the case cohort.

* **variance** - estimation of the influence-based variance, that follows the complete variance decomposition, for parameters such as log-relative hazard, cumulative baseline hazard or covariate specific pure-risk. Works with design weights or calibrated weights.

* **variance.missingdata** - estimation of the influence-based variance, that follows the complete variance decomposition, for parameters such as log-relative hazard, cumulative baseline hazard or covariate specific pure-risk, when covariate data is missing for individuals in the case cohort.

* **auxiliary.construction** - construction of the auxiliary variables proposed by Breslow et al. (Stat. Biosci., 2009) and Breslow et al. (IMS, 2013). 

* **calibration** - calibration of the design weights using the raking procedure. The Newton Raphson method is used for the optimization part.

* **product.covar.weight** - construction of the product of joint design weights and joint sampling indicators covariances, needed for the phase-two component of the variance (with design or calibrated weights).

* **cohort.generation** - generation of a cohort following the design used in the simulation study in Section 5 and Web Appendix E.

* **conf.interval** -construction of 95% confidence interval with parameter estimate and its variance estimate, assuming normality.

* **densityX1X3** and **densityX1** - used for the computation of the time-fixed baseline hazard in the simulations in Section 5 and Web Appendix E.

* **E.RH** and **E.RH.w** - computation of $\text{E} \lbrace \text{exp}\(\beta_1 X_1 + \beta_2 X_2 + \beta_3 X_3 \)\rbrace$ and $\text{E} \lbrace\text{exp}\(\beta_1 X_1 + \beta_2 X_2 + \beta_3 X_3 | \boldsymbol{W} = \boldsymbol{w} \)\rbrace$ , for $\boldsymbol{w}$ = {0, 1, 2, 3},  needed to compute the time-fixed baseline hazard and the sampling fractions to be used for an un-stratified and stratified sampling of the case-cohort in the simulations in Section 5 and Web Appendix E.

* **scenarios** - creation of a dataframe with the different scenarios to be investigated in the simulations in Section 5 and Web Appendix E.

* **scenarios.missingdata** - creation of a dataframe with the different scenarios to be investigated in the simulations in Web Appendix I (when covariate data is missing for individuals in the case cohort).
* **X.generation** - generation of a gaussian variable, with unit variance and given mean.
 
