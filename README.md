
<!-- README.md is generated from README.Rmd. Please edit that file -->
signatselect: Identifying signatures of selection
=================================================

[![Travis build status](https://travis-ci.org/benmarwick/signatselect.svg?branch=master)](https://travis-ci.org/benmarwick/signatselect) [![Binder](http://mybinder.org/badge_logo.svg)](http://beta.mybinder.org/v2/gh/benmarwick/signatselect/master?urlpath=rstudio)

The goal of signatselect is to provide two core functions useful for investigating change over time in artefact assemblages (and genetic time-series data):

-   `fit()` the frequency increment test as simple statistical test to aid in the detection and quantification of selective processes in the archaeological record. This is adapted directly from Feder, A. F., Kryazhimskiy, S., & Plotkin, J. B. (2014). Identifying signatures of selection in genetic time series. *Genetics*, 196(2), 509-522. <https://doi.org/10.1534/genetics.113.158220>.

-   `tsinfer()` to estimate the population size and the selection coefficient favoring one variant over another from time-series variant-frequency data. This is adapted from <https://github.com/bacovcin/tsinfer-R>

Installation
------------

You can install the development version of signatselect from GitHub with:

``` r
# install.packages("pak")
pak::pkg_install("benmarwick/signatselect")
```

Examples
--------

### The Frequency Increment Test

Here is an example of the `fit()`, the frequency increment test:

``` r
  # data slightly modified from Feder et al. Table 2 
  time <- c(415 , 505 , 585 , 665 , 745 , 825 , 910)
  freq <- c(0.06956522, 0.23125000, 0.62352941, 0.78494624, 0.93333333, 0.97979798, 0.98979592)
```

Let's take a look:

``` r
  plot(time, freq, type = 'b')
```

<img src="man/figures/README-fit-plot-1.png" width="50%" />

There's a trend of increasing frequencies, but is it a result of selection? Let's see:

``` r
library(signatselect)

    fit(
      time = time,
      v = freq
    )
#>   fit_stat      fit_p
#> 1 3.262457 0.02238466
```

The result of the FIT, with the low p-value, indicates that selection is occuring in this time series.

### Infer population size and selection coefficient from time-series variant-frequency data

Here is an example of `tsinfer()` to estimate the population size and the selection coefficient favoring one variant over another from time-series variant-frequency data. Here's some sample data: `tvec` contains sample times, `bvec` contains the number of samples of the focal variant (must be integers) and `nvec` containes total number of samples at each time point (must be integers).

``` r
  library(signatselect)
  # adapted from https://github.com/skryazhi/tsinfer
  tvec = c(0, 10, 20)
  bvec = c(2000, 4000, 6000)
  nvec = c(10000, 10000, 10000)
```

Now we compute the test result:

``` r
  tsinfer_output <- 
  tsinfer(
    tvec = tvec,
    bvec = bvec,
    nvec = nvec,
    verbose = FALSE
  )

# and take a look at the result 
tsinfer_output
#> $s.0
#> [1] 0
#> 
#> $alpha.0
#> [1] 80.02608
#> 
#> $f0.0
#> [1] 0.2
#> 
#> $LL.0
#> [1] 18.04008
#> 
#> $s
#> [1] 0.08925489
#> 
#> $alpha
#> [1] 16436.47
#> 
#> $f0
#> [1] 0.3659511
#> 
#> $LL
#> [1] 13.408
```

The selection coefficient for non-neutral model is in `s`, and so the value here is 0.0892549. The population size for non-neutral model is in `alpha`, and here is 1.643647410^{4}

An archaeological application
-----------------------------

Here is an example of using the FIT to identify pottery types that indicate selection. We are using frequencies of different decorative motifs in the Merzbach assemblage, Neolithic Germany (Crema et al. 2016; many other papers). We can load the data from the `evoarchdata` package on GitHub:

``` r
# pak::pkg_install("benmarwick/evoarchdata")
library(evoarchdata)
data("ceramics_lbk_merzbach")
```

Here's an overview of how each motif changes over time in this assemblage:

``` r
# get ordered factor of decoration types so we can order the plots nicely 
suppressPackageStartupMessages(library(tidyverse))

decoration_types <- 
names(ceramics_lbk_merzbach)[-1] %>%
  enframe() %>% 
  separate(value, into = c('a', 'b'), 2) %>% 
  mutate(b = parse_number(b)) %>% 
  arrange(b) %>% 
  unite(decorations, c(a,b), sep = "") %>% 
  pull(decorations)

# see how the freqs of each change over time
ceramics_lbk_merzbach_long <-
  ceramics_lbk_merzbach %>%
  gather(variable, value, -Phase) %>% 
  mutate(Phase = fct_relevel(Phase, ceramics_lbk_merzbach$Phase)) %>% 
  mutate(variable = fct_relevel(variable, decoration_types))

# plot
ggplot(ceramics_lbk_merzbach_long,
       aes(Phase,
           value)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  facet_wrap(~variable,
             scales = "free_y") +
  theme_minimal(base_size = 8)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

To prepare the data for the FIT we need to reshape it to a long form, compute frequency as ratio of count of type of a interest to all other types, and drop decoration types with less than threee time points (we need a min of three time points to compute the FIT).

``` r
# reshape data for each decoration type to go into test:
df <- ceramics_lbk_merzbach[ , 2:ncol(ceramics_lbk_merzbach)]
time <- utils:::.roman2numeric(ceramics_lbk_merzbach$Phase)

# compute frequency as ratio of count of type of interest to all other types
list_of_dfs <- vector("list", ncol(df))
names(list_of_dfs) <- names(df)

for(i in 1:ncol(df)){
  tmp <-
    data.frame(time = time,
               count_this_one = df[[i]],
               count_others = rowSums(df[, (seq_len(ncol(df)))[-i]   ]))
  
  # compute frequency
  tmp$frequency = with(tmp, count_this_one / count_others)
  
  # collect results and exclude rows with zero counts for this type i
  list_of_dfs[[i]] <- tmp[which(tmp$count_this_one != 0 ), ]
}

# we need a min of three time points to compute the FIT, so drop decoration types with less than 3
list_of_dfs_three_or_more <- 
  keep(list_of_dfs, ~nrow(.x) >= 3)
```

We can prepare safe version of the FIT so we can use it in loops without breaking out of the loop when there is an error for one iteration

``` r
fit_safely <- 
  safely(fit, 
         otherwise = data.frame(fit_stat = NA,
                                fit_p = NA))
```

Now we can compute the FIT for each pottery decoration type:

``` r
# apply test to each decoration type
df_fit_test_results <-
  list_of_dfs_three_or_more %>%
  bind_rows(.id = "type") %>%
  nest(-type) %>%
  mutate(fit_test = map(data,
                        ~fit_safely(time = .x$time,
                                    v =    .x$frequency))) %>%
  mutate(fit_p = map(fit_test, ~.x$result %>% bind_rows)) %>%
  unnest(fit_p) %>%
  mutate(sig = ifelse(fit_p <= 0.05, "selection", "neutral"))

ceramics_lbk_merzbach_long_sig <-
  ceramics_lbk_merzbach_long %>%
  left_join(df_fit_test_results %>% 
              select(type, sig), by = c("variable" = "type")) %>%
  mutate(Phase_num = utils:::.roman2numeric(as.character(Phase))) %>% 
  mutate(variable = fct_relevel(variable, decoration_types)) %>% 
  arrange(variable, Phase_num)

# plot to indicate which styles show selection and which do not. 
ggplot(ceramics_lbk_merzbach_long_sig,
       aes(Phase_num,
           value,
           colour = sig,
           shape = sig,
           group = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap(~variable,
             scales = "free_y") +
  theme_minimal(base_size = 8)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

We can subset each time series to see if we can identify episodes of selection among decoration types that might not show overall selection. First we need a function to isolate a data frame of rolling groups of three time points:

``` r
# Function to get df with rolling groups of 3
# from https://stackoverflow.com/a/5543632/1036500
df_with_rolling_idx <- function(df, window = 3) {
  nr <- nrow(df)
  w <- window       # window size
  i <- 1:nr         # indices of the rows
  iw <-
    embed(i, w)[, w:1]   # matrix of rolling-window indices of length w
  wnum <- rep(1:nrow(iw), each = w)   # window number
  inds <-
    i[c(t(iw))]         # the indices flattened, to use below
  zw <- sapply(df, '[', inds)
  zw <- transform(data.frame(zw), w = wnum)
  return(zw)
}
```

Now we can compute the FIT on sections of the time series for each decoration type to identify time coordinates where selection has occurred, even when the overall series does not indicate selection:

``` r
n <- 5
merzbach_long_sig_mid_time_point <- 
list_of_dfs %>% 
  keep(., ~nrow(.x) >= n) %>% # need at least 4 rows of a type
  bind_rows(.id = "type") %>% # to get rolling window of 3
  nest(-type) %>%
  mutate(rolled = map(data, df_with_rolling_idx)) %>% 
  unnest(rolled) %>% 
  mutate(unid = str_glue('{type}_{w}')) %>% 
  nest(-unid) %>% 
  mutate(fit_test = map(data,
                        ~fit_safely(time = .x$time,
                                    v =    .x$frequency))) %>%
  mutate(fit_p = map(fit_test, ~.x$result %>% bind_rows)) %>%
  unnest(fit_p) %>%
  mutate(sig = ifelse(fit_p <= 0.05, "selection", "neutral")) %>% 
  unnest(data) %>%                            
  group_by(unid) %>%                          
  slice(ceiling(n()/2)) %>% 
  right_join(list_of_dfs %>% 
               keep(., ~nrow(.x) >= n) %>% # need at least 4 rows of a type
               bind_rows(.id = "type")) 
#> Joining, by = c("type", "time", "count_this_one", "count_others", "frequency")

# make type a factor so we can order the plots nicely 
merzbach_long_sig_mid_time_point$type <- 
  fct_relevel(merzbach_long_sig_mid_time_point$type,
              decoration_types[decoration_types %in% merzbach_long_sig_mid_time_point$type])

# plot with overall time-series results also
# harmonize some variable names first

ceramics_lbk_merzbach_long_sig_to_plot_with_others <- 
  ceramics_lbk_merzbach_long_sig %>% 
  rename( time = Phase_num,
          count_this_one = value,
          type = variable) %>% 
  filter(count_this_one != 0) %>% 
  arrange(type, time) %>% 
  mutate(type = fct_relevel(type, decoration_types))

# here we have the plot showing overall selection, and point-wise selection
ggplot()  +
  geom_line(data = merzbach_long_sig_mid_time_point %>% 
              filter(sig == "selection"),
            aes(time,
                count_this_one,
                group = type),
            size = 5,
            colour = "grey80",
            lineend = "round") +
  geom_point(data = merzbach_long_sig_mid_time_point %>% 
               filter(sig == "selection") %>% 
               group_by(type)  %>% 
               filter(n() > 2),
             aes(time,
                 count_this_one,
                 group = type),
             size = 5,
             colour = "grey80") +
  geom_line(data = ceramics_lbk_merzbach_long_sig_to_plot_with_others,
            aes(time,
                count_this_one, 
                group = type,
                colour = sig)) +
  geom_point(data = ceramics_lbk_merzbach_long_sig_to_plot_with_others,
            aes(time,
                count_this_one, 
                group = type,
                colour = sig,
                shape = sig))  +
  facet_wrap( ~ type, scales = "free_y") +
  theme_minimal(base_size = 8)
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to
#> adjust the group aesthetic?
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Please note that the `signatselect` project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.
