# drconfseq: Doubly robust confidence sequences

This package implements doubly robust confidence sequences which provide
safe, anytime-valid inference for the average treatment effect (ATE)
from sequential data. These allow you to

1.  Update inferences about the ATE as new data become available;
2.  Continuously monitor studies and adaptively stop or continue them
    for any reason; and
3.  Control the type-I error at all stopping times, including random or
    data-dependent times.

Furthermore, these confidence sequences can be used in both randomized
sequential experiments and observational studies with no unmeasured
confounding.

The main methodological reference is

-   Waudby-Smith, I., Arbour, D., Sinha, R., Kennedy, E. H., &
    Ramdas, A. (2021). Doubly robust confidence sequences for sequential
    causal inference. arXiv preprint
    [arXiv:2103.06476](https://arxiv.org/pdf/2103.06476.pdf).

## Installation

You can install the package using `devtools`:

``` r
# Install devtools if you do not have it:
install.packages("devtools")

# Install the package using devtools::install_github 
devtools::install_github("wannabesmith/drconfseq")
```

## Example usage

Let’s jump into a simple example taken from [“Doubly-robust confidence
sequences for sequential causal
inference”](https://arxiv.org/pdf/2103.06476.pdf). Suppose that we wish
to estimate the ATE in an observational setting with no unmeasured
confounding.

``` r
library(drconfseq)

# Optional, can be used to speed up computations
library(parallel)
```

First, we will generate *n* = 10<sup>4</sup> observations each with 3
real-valued covariates from a trivariate Gaussian.

``` r
n <- 10000
d <- 3
X <- cbind(1, matrix(rnorm(n*d), nrow = n))
```

Randomly assign subjects to treatment or control groups with equal
probability.

``` r
treatment <- rbinom(n = n, size = 1, prob = 1 / 2)
```

Define the regression function, and the target parameter (which we will
ensure is the average treatment effect by design),

*A**T**E* := 1.

``` r
beta_mu <- c(1, -1, -2, 3)

mu <- function(x)
{
  beta_mu %*% c(1, x[1]^2, sin(x[2]), abs(x[3]))
}

ATE <- 1
```

Finally, generate outcomes *y*<sub>1</sub>, …, *y*<sub>*n*</sub> as

*y*<sub>*i*</sub> := *f*<sup>⋆</sup>(*x*<sub>*i*, 1</sub>, *x*<sub>*i*, 2</sub>, *x*<sub>*i*, 3</sub>) + *A**T**E* ⋅ *t**r**e**a**t**m**e**n**t*<sub>*i*</sub> + *ϵ*<sub>*i*</sub>,

where *ϵ*<sub>*i*</sub> ∼ *t*<sub>5</sub> are drawn from a
*t*-distribution with 5 degrees of freedom.

``` r
y <- apply(X, MARGIN=1, FUN=mu) + treatment*ATE + rt(n, df=5)
```

We will build three estimators for the ATE with increasing degrees of
complexity.

1.  **Unadjusted**: uses the constant function 0 for the outcome
    regression. This is equivalent to an inverse-probability-weighted
    estimator.
2.  **Parametric**: estimates the outcome regression with a linear
    model.
3.  **Super Learner**: estimates the outcome regression with a “Super
    Learner” weighted average of linear models and other machine
    learning algorithms including regression splines, random forests,
    and generalized additive models.

Since we are in a randomized experiments, all three estimators will
consistently estimate and provide valid sequential inference for the
ATE. However, since the underlying regression function is nonlinear, we
expect the Super Learner to outperform the Parametric to outperform the
Unadjusted.

``` r
# Get SuperLearner prediction function for $\mu^1$.
sl_reg_1 <- get_SL_fn(SL.library = c("SL.earth", "SL.gam", "SL.glm", "SL.ranger"))

# Do the same for $\mu^0$.
sl_reg_0 <- sl_reg_1 

# Get Parametric regression function for $\mu^1$.
parametric_reg_1 = get_SL_fn(SL.library = "SL.glm")

# Do the same for $\mu^0$.
parametric_reg_0 = get_SL_fn(SL.library = "SL.glm")

# Compute the confidence sequence at logarithmically-spaced time points
times <- unique(round(logseq(500, 10000, n = 10)))
alpha <- 0.1

# Set n_cores to 1 if you are not using the parallel package
n_cores <- detectCores()

confseq_SL <- confseq_ate(y, X, treatment, regression_fn_1 = sl_reg_1,
                          regression_fn_0 = sl_reg_0,
                          propensity_score_fn = function(y, X, newX){1/2},
                          t_opt = 1000, alpha=alpha,
                          times=times, n_cores = n_cores, cross_fit = TRUE)
confseq_glm <- confseq_ate(y, X, treatment, regression_fn_1 = parametric_reg_1,
                           regression_fn_0 = parametric_reg_0,
                           propensity_score_fn = function(y, X, newX){1/2},
                           t_opt = 1000, alpha=alpha,
                           times=times, n_cores = n_cores, cross_fit = TRUE)

confseq_unadj <- confseq_ate_unadjusted(y = y, treatment = treatment,
                                        propensity_score = 1/2,
                                        t_opt = 1000, alpha = alpha,
                                        times = times)
```

Let us now plot the confidence sequence.

``` r
library(ggplot2)

plot_df <- data.frame(t = rep(times, 6), 
                      bound = c(confseq_SL$l, confseq_SL$u,
                                confseq_glm$l, confseq_glm$u,
                                confseq_unadj$l, confseq_unadj$u),
                      upper_lower = c(rep('l', length(confseq_SL$l)),
                                      rep('u', length(confseq_SL$u)),
                                      rep('l', length(confseq_glm$l)),
                                      rep('u', length(confseq_glm$u)),
                                      rep('l', length(confseq_unadj$u)),
                                      rep('u', length(confseq_unadj$u))),
                      Model = c(rep('Super Learner', 2*length(confseq_SL$l)),
                               rep('Parametric', 2*length(confseq_glm$l)),
                               rep('Unadjusted', 2*length(confseq_unadj$l))))

plt_cs <- 
  ggplot() + 
  geom_line(aes(x=t, y=bound, color=Model, linetype=Model), 
            data=plot_df[plot_df$upper_lower=='l',]) +
  geom_line(aes(x=t, y=bound, color=Model, linetype=Model), 
            data=plot_df[plot_df$upper_lower=='u',]) +
  geom_hline(yintercept=ATE, linetype='dashed', color='gray') +
  ylab('Confidence sequences for the\naverage treatment effect') + xlab('Time') +
  scale_x_log10(breaks = scales::trans_breaks("log10",
                                              function(x) 10^x),
                labels = scales::trans_format("log10",
                                              scales::math_format(10^.x))) +
  annotation_logticks(colour = "grey", side = "b") + 
  theme_minimal() + 
  theme(text = element_text(family="serif")) +
  scale_color_brewer(palette="Set2")

plt_cs
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Run unit tests

First, clone the repository:

``` zsh
git clone git@github.com:WannabeSmith/drconfseq.git
```

Then in an R console, run

``` r
devtools::check("path/to/drconfseq")
```

## Reproduce the paper’s plots

See the README in [this folder](paper_plots).

## Credits

The methods implemented in this package were developed by [Ian
Waudby-Smith](https://ian.waudbysmith.com), [David
Arbour](https://darbour.github.io/), [Ritwik
Sinha](https://research.adobe.com/person/ritwik-sinha/), [Edward H.
Kennedy](http://ehkennedy.com), and [Aaditya
Ramdas](http://stat.cmu.edu/~aramdas).
