# TSINFER

Infer population size and selection coefficient from time-series allele-frequency data


## Requirements
  * OCaml langugage compiler: http://ocaml.org/
  * GNU Scientific Library (GSL): http://www.gnu.org/software/gsl/
  * GSL Interface for OCaml language: http://oandrieu.nerim.net/ocaml/


## Installation
Untar the archive into the desired directory::
```
  tar -zxf tsinfer.tar.gz
```

Compile required modules and the executable:
```
  make all
```


## Execution
Execution syntax is simple
```
  ./tsinfer infile [OPTIONS]
```
infile -- file containing the time series data. The file should contain three lines, each tab separated, with equal number of entries. The first line contains sample times, the second line contains the number of samples of the focal allele (must be integers) and the third line containes total number of samples at each time point (must be integers). For example

0  10      20  2000    4000    6000
10000   10000   10000

  OPTIONS:
  -f -- treat sampled frequencies as actual frequencies
  -func_eval:XX -- number of function evalutations in the Monte Carlo integration (posiitve integer). Default : 3e+04
  -h -- help
  -maxiter:XX -- maximum number of iterations (positive integer). Default : 200
  -maxterms:XX -- maximum number of terms in the series (positive integer). Default 2000
  -neut -- evaluate only the likelihood of the data under neutrality
  -prec:XX -- precision (for the size of the simples). Default : 1.0e-06
  -start:s,alpha -- starting guess for s (float) and alpha (positive float). Default : (0, 1.000e+03)
  -stds:XX -- number of standard deviations around the mean of the beta distribution for the approximation of the intergral (positive integer). Default : 3.0
  -step:s,alpha -- step size for s, and alpha (all positive floats). Default : (0.005, 5.000e+01)
  -plot:filename -- instead of finding the maximum likelihood value, plot the likelihood function profile at points specified in the file 'filename'. Each line of the file must consist of 2 tab-separated numbers corresponding to the coordinates (s, alpha). See examples. The 'start', and 'step' options are ignored.
  -v -- verbose mode

The output of the execution is a list of 10 values, e.g.,   STATUS = SUCCESS
MLS = 0.08927
MLALPHA = 1.635e+04
MNLOGL = 20.21
AIC2 = 44.419
S0ALPHA = 7.999e+01
S0NLOGL = 24.837
AIC1 = 51.675
LLR = 9.2554
CHI2P = 2.4e-03

STATUS indicates whether the program successfully found the ML values of selection coefficient s and the population size alpha. If procedure is not successful the failure mode is indicated.
MLS is the ML value of the selection coefficient s in the two-parameter model
MLALPHA is the ML value of the population size alpha in the two-parameter model
MNLOGL is the obtained minimum value of the negative log-likelihood in the two-parameter model
AIC2 is the value of the Akaike information criterion for the two-parameter model
S0ALPHA is the ML value of population size alpha in the one-parameter model (i.e., when s = 0).
S0NLOGL is the obtained minimum value of the negative log-likelihood in the one-parameter model
AIC1 is the value of the Akaike information criterion for the one-parameter model
LLR is the twice log likelihood ratio statistic, i.e., 2 * (S0NLOGL - MNLOGL)
CHI2P is the right-tailed chi2 P-value for the likelihood ratio statistic

## Examples
**Example 1.** Obtain the maximum likelihood parameters for time series in file test.in, taking sampling into account
```
  ./tsinfer test.in
```
Note that, by default, the procedure of calculating the likelihood of data takes binomial sampling into account. However, this procedure was not extensively tested and may contain errors. It is recommende to use instead the procedure that treats sampled frequecies as exact. See next example below.


**Example 2.** Obtain the maximum likelihood parameters for time series in file test.in, treating sampled frequencies as exact. Note that the input file still needs to have three lines with second and third lines containing integer values.
```
  ./tsinfer test.in -f
```

**Example 3.** Obtain the values of the negative log-likelihood function for time series in file test.in at parameter values specified in file grid, treating sampled frequencies as exact
```
  ./tsinfer test.in -f -plot:grid
```
