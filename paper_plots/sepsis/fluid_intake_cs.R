library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)

library(drconfseq)

source("clean_data.R")

sepsis <- fread("./data/sepsis_patients.csv",
  na.strings = "NULL"
) %>%
  clean_sepsis_data()

start_time <- 1000
times <- unique(round(logseq(start_time, nrow(sepsis), n = 30)))
t_opt <- 1500
alpha <- 0.1

cs_fluid_intake <- asymptotic_confseq(x = sepsis$fluids_24h_l, t_opt = t_opt, alpha = alpha)

cs_fluid_intake$l <- cs_fluid_intake$l[times]
cs_fluid_intake$u <- cs_fluid_intake$u[times]
muhat_fluid_intake <- cumul_mean(sepsis$fluids_24h_l)[times]

r_data_dir <- "./"
save(cs_fluid_intake,
  muhat_fluid_intake,
  start_time,
  times,
  file = paste(r_data_dir, "fluid_intake_cs.RData", sep = "")
)

## save(cs_fluid_intake, file=paste(r_data_dir, "fluid_intake_cs.RData", sep=""))
