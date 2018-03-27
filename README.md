# NBARX #

This repository contain all the R-code used for the master thesis "Maximum likelihood estimation of a negative binomial autoregression with exogenous co-variates (NBARX) - applied to Danish foreclosures".

The script `func.R` is a standalone function library containing likelihood functions, filters and more for estimating Poisson Autoregression (PAR), Negative Binomial Autoregression (NBAR), Poisson Autoregression with exogenous variables (PARX) and Negative Binomial Autoregression with exogenous varables (NBARX).

The purely autoregressive models PAR and NBAR are fitted to the danish monthly series of,

- foreclosures
- bankrupties
- marriages
- divorces

both seasonally adjusted and non-adjusted version, so eight series in total.

The two models containing exogenous variables, PARX and NBARX, model the monthly danish foreclosures counts
as a function of the exogenous variables,

- House price index of one family houses
- Net unemployment
- Divorce counts


### Dependencies
#### The function library func.R
`library(sfsmisc)`  
`library(numDeriv)`  
`library(zoo)`

#### The rest of the thesis code
`library(readxl)`  
`library(zoo)`  
`library(RColorBrewer)`  
`library(season)`

### The different scripts
#### func.R
This is a standalone function library which contains mainly three type of functions
intended for estimation and simulation of linear generalized autoreggresive integer models (INGAR), the
four specific models are PAR, NBAR, PARX and NBARX models.
The three functions type are:

- simulation functions: `sim.PAR`, `sim.NBAR`, `sim.PARX` and `sim.NBARX`
- likelihood functions: `loglike.PAR`, `loglike.NBAR`, `loglike.PARX` and `loglike.NBARX`
- one step forecast filter: : `filter.PAR`, `filter.NBAR`, `filter.PARX` and `filter.NBARX`

The function library also contains a function for calculating the theoretical autocovariance and 
autocorrelation functions (ACF), based on the ARMA(X)-representation of the models. The `acf`-functions
can also be used to calculate ACF's of GARCH and GARCH-X models.

#### est y.R and est X.R
These scripts make all the estimation on the 8 data series fetched from Statistics Denmark, and saves the
estimated parameters, the filtered series and some test results in RDS-format (R-data single) 
in the `est_output` folder.

`est y.R` fits the PAR and NBAR model on all 8 series.

`est X.R` fits the PARX and NBARX models to foreclosure counts, with changes in house prices, unemployment and divorces
as explanatory variables.

#### plot y.R and plot X.R
Makes all the plots for the thesis, based on the estimation output in `est_output` and saves the plots in pdf format 
in the folder `plots`.

#### get_data_StatsBank.R
Downloads data through Statistics Denamarks API, 
[www.dst.dk/en/Statistik/statistikbanken/api](http://www.dst.dk/en/Statistik/statistikbanken/api),
by posting the JSON-strings in the folder POST_json, through the API. The data is then
pre-processed (e.g. som series are stitched together to get more observations), for some series
seasonality adjustment is made with R-package `season` and finally the data is saved in the folder `data`.


