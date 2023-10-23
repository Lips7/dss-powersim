


# (PART) Appendix {-}

# Single Run Benchmark

In this section, we will measure the performance differences among r, python and stata, specifically the `single_run` function of each software,
which is the atomic function of whole power simulation program.
But, due to the differences of parameters selection, the measurement is not technically subjective, the benchmark results here are for reference only.

## R

In implementation with [R](./r.html), we used all parameters ($\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$) to generate simulated data
and used the function `lmer()` from the package `{lmerTest}` to perform linear mixed effects analyze.

Here are the definitions of functions [sim_data](./r.html#data-simulation-automated) and [single_run](./r.html#power-calculation-single-run).



We'll use `benchmark` function from package `{rbenchmark}` to measure the running time of function `single_run`.


```r
library(rbenchmark)

result <- within(
  benchmark(
    single_run = single_run(),
    replications = 100,
    columns = c("test", "elapsed", "replications")
  ),
  {
    average <- elapsed / replications
  }
)

r_average <- result[1, "average"]

sprintf("Average running time is %f s", r_average)
```

```
## [1] "Average running time is 0.232100 s"
```

## Python

In implementation with [Python](./python.html), we used all parameters ($\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$) to generate simulated data
and used the function `mixedlm()` from package `{statsmodels}` to perform linear mixed effects analyze.

Here are the definitions of functions [sim_data](./python.html#data-simulation-automated) and [single_run](./python.html#power-calculation-single-run).





We'll use `timeit` function from the package `{timeit}` to measure the running time of function `single_run`.


```python
import timeit

result = timeit.timeit(single_run, number=100)
python_average = result/100

print(f"Average running time is {python_average} s")
```

```
## Average running time is 1.1418413130000409 s
```

## Stata

In implementation with [Stata](./stata.html), we only used 4 parameters ($\beta_0$, $\beta_1$, $\tau_0$, $\sigma$) to generate simulated data
and used the function `mixed` to perform linear mixed effects analyze.

Here are the definitions of functions [sim_data](./stata.html#data-simulation-automated) and [single_run](./stata.html#power-calculation-single-run).





In stata, we have to implement the benchmark code by our own.


```stata
capture program drop tstart
program tstart, rclass
  timer clear 1
  timer on 1
end

capture program drop tend
program tend, rclass
  timer off 1
  qui timer list 1
  scalar r = r(t1)
end

tstart
quietly {
    forval i = 1/100 {
        single_run 25 15 15 60 5 7 8
    }
}
tend
local result = scalar(r)
local stata_average = r / 100
di "Average running time is " `stata_average' " s"

quietly {
    file open result using "./data/stata_average.txt", write replace
    set more off

    file write result "`stata_average'"
    file close result
    set more on
}
```

```
Average running time is .23617 s
```

## Benchmark Results

Finally, we get all the results of r, python and stata, let's show them with a table.

Table: (\#tab:unnamed-chunk-9) Benchmark results.


Software   Parameters                                                               Average Time 
---------  -----------------------------------------------------------------------  -------------
R          $\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$   0.23         
Python     $\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$   1.14         
Stata      $\beta_0$, $\beta_1$, $\tau_0$, $\sigma$                                 0.24         
