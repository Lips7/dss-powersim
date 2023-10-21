


# (PART) Appendix {-}

# Software Comparison

Each of part has software basic descriptions, model descriptions, and `single_run` function benchmark results (`single_run` is the atomic function of power simulation program).

Note: The benchmark results presented below are based on a laptop, here are the device details.
- Model: Microsoft Laptop 3
- Windows Edition: Windows 11 Home Insider Preview 25977.1000
- CPU: Intel(R) Core(TM) i7-1065G7 CPU @ 1.30GHz 4 cores
- RAM: 16 GB

Packages version
- R: 4.3.1
  - lmerTest: 3.1
- Python: 3.11.6
  - numpy: 1.26.1
  - pandas: 2.1.1
  - statsmodels: 0.14.0
- Stata: Stata17MP

## R

R is a dynamic and high-efficient programming language for statistical computing and graphics.
In implementation with [R](./r.html), we used all parameters ($\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$) to generate simulated data
and used the function `lmer()` from the package `{lmerTest}` to perform linear mixed effects analyze.



First we need to define the functions `sim_data` and `single_run`.


```r
# set up the custom data simulation function
sim_data <- function(
  n_subj     =  25,   # number of subjects
  n_pop      =  15,   # number of pop songs
  n_rock     =  15,   # number of rock songs
  beta_0     =  60,   # mean for pop genre
  beta_1     =   5,   # effect of genre
  omega_0    =   3,   # by-song random intercept sd
  tau_0      =   7,   # by-subject random intercept sd
  tau_1      =   4,   # by-subject random slope sd
  rho        =   0.2, # correlation between intercept and slope
  sigma      =   8    # residual (standard deviation)
) {
  # simulate a sample of songs
  songs <- tibble(
    song_id = seq_len(n_pop + n_rock),
    category = rep(c("pop", "rock"), c(n_pop, n_rock)),
    genre_i = rep(c(0, 1), c(n_pop, n_rock)),
    O_0i = rnorm(n = n_pop + n_rock, mean = 0, sd = omega_0)
  )

  # simulate a sample of subjects
  subjects <- faux::rnorm_multi(
    n = n_subj,
    mu = 0,
    sd = c(tau_0, tau_1),
    r = rho,
    varnames = c("T_0j", "T_1j")
  ) |>
  mutate(subj_id = seq_len(n_subj))

  # cross subject and song IDs
  crossing(subjects, songs) |>
    mutate(e_ij = rnorm(n(), mean = 0, sd = sigma),
         liking_ij = beta_0 + T_0j + O_0i + (beta_1 + T_1j) * genre_i + e_ij) |>
  select(subj_id, song_id, category, genre_i, liking_ij)
}

# set up the power function
single_run <- function(...) {
  # ... is a shortcut that forwards any additional arguments to sim_data()
  dat_sim <- sim_data(...)
  mod_sim <- suppressWarnings({ suppressMessages({ # suppress singularity messages
    lmerTest::lmer(liking_ij ~ 1 + genre_i + (1 | song_id) + (1 + genre_i | subj_id), data = dat_sim)
  })})
  broom.mixed::tidy(mod_sim)
}
```

Then we'll use `benchmark` from the package `{rbenchmark}` to measure the running time of function `single_run`.


```r
library(rbenchmark)

result <- within(
    benchmark(
        single_run = single_run(),
        replications = 100,
        columns = c('test', 'elapsed', 'replications')
    ),
    { average = elapsed/replications }
)

r_average <- result[1, 'average']

sprintf("Average running time is %f s", r_average)
```

```
## [1] "Average running time is 0.175900 s"
```

## Python

Python is a high-level, general-purpose programming language.
In implementation with [Python](./python.html), we used all parameters ($\beta_0$, $\beta_1$, $\omega_0$, $\tau_0$, $\tau_1$, $\rho$, $\sigma$) to generate simulated data
and used the function `mixedlm()` from the package `{statsmodels.formula.api}` to perform linear mixed effects analyze.





Like we did in R, we define functions `sim_data` and `single_run`.


```python
def sim_data(n_subj=25, n_pop=15, n_rock=15, beta_0=60, beta_1=5, omega_0=3, tau_0=7, tau_1=4, rho=0.2, sigma=8):
    songs = pd.DataFrame({
        'song_id': np.arange(n_pop + n_rock),
        'category': np.repeat(["pop", "rock"], [n_pop, n_rock]),
        'genre_i': np.repeat([0, 1], [n_pop, n_rock]),
        'O_0i': np.random.normal(0, omega_0, n_pop + n_rock)
    })

    random_effects = np.random.multivariate_normal([0, 0], [[tau_0**2, rho*tau_0*tau_1], [rho*tau_0*tau_1, tau_1**2]], n_subj)
    subjects = pd.DataFrame(random_effects, columns=['T_0j', 'T_1j'])
    subjects['subj_id'] = np.arange(1, n_subj + 1)

    trials = subjects.assign(key=1).merge(songs.assign(key=1), on='key').drop(columns='key')
    trials['e_ij'] = np.random.normal(0, sigma, len(trials))
    trials['liking_ij'] = beta_0 + trials['T_0j'] + trials['O_0i'] + (beta_1 + trials['T_1j']) * trials['genre_i'] + trials['e_ij']

    return trials[['subj_id', 'song_id', 'category', 'genre_i', 'liking_ij']]

def single_run(n_subj=25, n_pop=15, n_rock=15, beta_0=60, beta_1=5, omega_0=3, tau_0=7, tau_1=4, rho=0.2, sigma=8):
    dat_sim = sim_data(n_subj, n_pop, n_rock, beta_0, beta_1, omega_0, tau_0, tau_1, rho, sigma)

    dat_sim['groups'] = 1
    mod_sim = smf.mixedlm('liking_ij ~ 1 + genre_i', groups=dat_sim['groups'],
                          vc_formula={'song_id':'0 + C(song_id)', 'subj_id':'0 + C(subj_id)', 'genre_i': '0 + C(subj_id):genre_i'},
                          re_formula='0', data=dat_sim).fit()

    df = mod_sim.summary().tables[1]
    df['p_value'] = mod_sim.pvalues
    return df[['Coef.', 'Std.Err.', 'p_value']]
```

And then we'll use `timeit` function from the package `{timeit}` to measure the running time of function `single_run`.


```python
import timeit

result = timeit.timeit(single_run, number=100)
python_average = result/100

print(f"Average running time is {python_average} s")
```

```
## Average running time is 0.8212643380000009 s
```

## Stata

Stata is a general-purpose statistical software package developed by StataCorp for data manipulation, visualization, statistics, and automated reporting.
In implementation with [Stata](./stata.html), we only used 4 parameters ($\beta_0$, $\beta_1$, $\tau_0$, $\sigma$) to generate simulated data
and used the function `mixedl` to perform linear mixed effects analyze.





Like we did in R, we define functions `sim_data` and `single_run`.


```stata
capture program drop sim_data
program define sim_data
	args n_subj n_pop n_rock beta_0 beta_1 tau_0 sigma

	clear
	local n_all = `n_pop' + `n_rock'
	set obs `n_all'
	gen song_id = _n
	gen category = "pop"
	replace category = "rock" if song_id > `n_pop'
	gen genre_i = 0
	replace genre_i = 1 if song_id > `n_pop'
	gen key = 1
	save "./data/songs.dta", replace

	clear
	set obs `n_subj'

	gen t0 = rnormal(0, `tau_0')
	gen subj_id = _n
	gen key = 1
	save "./data/subjects.dta", replace

	use "./data/subjects.dta"
	cross using "./data/songs.dta"
	drop key
	sort subj_id song_id
	gen e_ij = rnormal(0, `sigma')

	gen liking_ij = `beta_0' + t0 + `beta_1' * genre_i + e_ij
	keep subj_id song_id category genre_i liking_ij
end

capture program drop single_run
program define single_run, rclass
	args n_subj n_pop n_rock beta_0 beta_1 tau_0 sigma

	clear
	sim_data `n_subj' `n_pop' `n_rock' `beta_0' `beta_1' `tau_0' `sigma'
	mixed liking_ij genre_i || subj_id:, noretable nofetable noheader nogroup

    estimates clear
	estimates store model_results

	matrix coefficients = e(b)
	matrix std_errors = e(V)
	matrix p_values = e(p)

	return scalar coef = coefficients[1, 1]
	return scalar std_err = std_errors[1, 1]
	return scalar p_value = p_values[1, 1]
end
```

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
Average running time is .23042 s
```

## Results

Finally, we get all the results of r, python and stata, let's show them with a table. And the conclusion is r > stata >> python.

Table: (\#tab:unnamed-chunk-12) Benchmark results.


lang     parameters                                 average 
-------  -----------------------------------------  --------
r        All parameters                             0.18    
python   All parameters                             0.82    
stata    $\beta_0$, $\beta_1$, $\tau_0$, $\sigma$   0.23    
