# TSINFER

An R package implementation of code to infer population size and selection coefficient from time-series allele-frequency data


## Requirements
  * R statistical language
  * nloptr package for R


## Installation
 In R:

 install.package('devtools')
 install_github('bacovcin/tsinfer-R',quick=T)

## Execution
 Use the tsinfer function.

 Essential arguments are tvec (time point labels for time series starting at 0), bvec (number of new alleles at each time point), and nvec (total number of samples at teach time point).

### Output
The output of the execution is the following list:

s = selection coefficient for non-neutral model
alpha = population size for non-neutral model
f0 = initial frequency for the logistic in non-neutral model
LL = log-likelihood of non-neutral model
s.0 = 0 (selection coefficient for the neutral model)
alpha.0 = population size for the neutral model
f0.0 = initial frequency for the logistic in neutral model
LL.0 = log-likelihood of neutral model
