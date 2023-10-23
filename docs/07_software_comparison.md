


# (PART) Appendix {-}

# Software Comparison

## Summary

In [R](./r.html), we implemented the power simulation process completely, we used all parameters
($\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$) to generate simulated data
(see functions [sim_data](./r.html#data-simulation-automated) and [single_run](./r.html#power-calculation-single-run)).

In [Python](./python.html), we also used all parameters ($\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$) to generate simulated data
(see functions [sim_data](./python.html#data-simulation-automated) and [single_run](./python.html#power-calculation-single-run)).
But there are two slight differences between r and python:

- Simulation of the sampling of subjects:
    - R: we used the function `rnorm_multi` form package `{faux}` to generate a table of `n` simulated value from a multivariate normal distribution by specifying
    the mean ($\mu = 0$) and **standard deviations** ($sd = (\tau_0, \tau_1)$), plus the **correlations** ($r = \rho$).
    - Python: we used the function `multivariate_normal` from package `{numpy}` to generate a table of `n` simulated values from a multivariate normal distribution by specifying
    the means ($\mu = [0, 0]$) and **covariance matrix** ($cov = [[\tau_0^2, \rho * \tau_0 * \tau_1], [\rho * \tau_0 * \tau_1, \tau_1^2]]$).

- Mixed Effects Models:
    - R: we indicated the correlation between the random intercept and random slope of subjects in the formula of `lmer` function `liking_ij ~ 1 + genre_i + (1 | song_id) + (1 + genre_i | subj_id)`.
    - Python: due to the inability of the function `mixedlm` from package `{statsmodels}`, we didn't indicate the correlations in the model.

In [Stata](./stata.html), we only used 4 parameters ($\beta_0$, $\beta_1$, $\tau_0$, $\sigma$) to generate simulated data.
So there are noticeable changes in the code of data simulation and mixed effects analyze process
(see functions [sim_data](./stata.html#data-simulation-automated) and [single_run](./stata.html#power-calculation-single-run)).

## Absolute Running Time

Note: All the results presented below are based on a laptop, here are the device details.

Table: (\#tab:unnamed-chunk-1) Basic Information

                                          Info                        
-----------------  ---------------------------------------------------
      Model                        Microsoft Laptop 3                 
 Windows Edition       Windows 11 Home Insider Preview 25977.1000     
       CPU          Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz 4 cores 
       RAM                              DDR4 16GB                     

Now we will list each absolute running time of each software, there are a lot of differences among them, here are the details.

Table: (\#tab:unnamed-chunk-2) Absolute Running Time


 Software                                                                              Info                                                                                                                                   Version                                                                                       Parameters                                  Replicates    Absolute Time 
----------  -----------------------------------------------------------------------------------------------------------------------------------------------------------  ------------------------------------------------------------------------------------------------------------------  ------------------------------------------------------------------------  ------------  ---------------
    R                                     R is a dynamic and high-efficient programming language for statistical computing and graphics.                                                      $\substack{\text{R: 4.3.1} \\ \text{lmerTest: 3.1} \\ \text{lme4: 1.1}}$                        $\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$       100            4m13s     
  Python                                                   Python is a high-level, general-purpose programming language.                                                  $\substack{\text{Python: 3.11.6} \\ \text{numpy: 1.26.1} \\ \text{pandas: 2.1.1} \\ \text{statsmodels: 0.14.0}}$    $\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$        30           21m27s     
  Stata      Stata is a general-purpose statistical software package developed by StataCorp for data manipulation, visualization, statistics, and automated reporting.                                          $\substack{\text{Stata: Stata17MP}}$                                                         $\beta_0$, $\beta_1$, $\tau_0$, $\sigma$                       30           6min7s     
