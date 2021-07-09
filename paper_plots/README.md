# Reproduce paper plots

Here we describe how to reproduce the plots in [Doubly robust confidence sequences for sequential causal inference](https://arxiv.org/pdf/2103.06476.pdf). This is broken up into three steps:
1. Run simulations,
2. [Analyze sepsis data](todo) (optional), and
3. Generate plots from the output of steps 1 and 2.

The sepsis data analysis is labeled as "optional" because it requires access to, and extracting from, the [MIMIC-III](https://mimic.mit.edu/docs/iii/) database, which requires significant memory and computational resources. For instructions on how to extract data from MIMIC-III, see the README found in [this folder](sepsis).

## 1. Run simulations

First, make sure that `sequential.causal` is installed in your R environment:

```R
# Install devtools from CRAN if you don't have it: install.packages("devtools")
devtools::install_github("wannabesmith/sequential.causal")
```

Now, clone the git repository:

```zsh
git clone git@github.com:WannabeSmith/sequential.causal.git
```

Navigate to the simulations folder:

```zsh 
cd sequential.causal/paper_plots/simulations/
```

Run all simulations using the shell script. For example,

```zsh
# Use the shell of your choice, e.g. zsh.
# Note, this could take a while...
zsh generate_plot_RData.sh
```

Alternatively, you can run each simulation separately. For example,
```zsh
cd observational_study_SL/
Rscript observational_study_SL.R
```

and so on.

## 2. Analyze sepsis data (optional)

See the README in [this folder](sepsis).

## 3. Generate plots

Open `sequential.causal/paper_plots/plots_for_paper.Rmd` with RStudio and run the chunks sequentially. Alternatively you can simply knit the Rmd file.

Note that if you did not complete step 2, the final chunk will not run correctly (and neither will a full knit).
