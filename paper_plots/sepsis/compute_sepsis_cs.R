library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)

library(drconfseq)

clean_sepsis_data <- function(data) {
  data %>%
    lazy_dt(immutable = TRUE) %>%
    mutate(
      # Simplifying ethnicity column according to NIH guidelines
      # https://grants.nih.gov/grants/guide/notice-files/not-od-15-089.html
      ethnicity_NIH =
        case_when(
          grepl("AMERICAN INDIAN", ethnicity) ~ "native",
          grepl("ASIAN|MIDDLE EASTERN", ethnicity) ~ "asian",
          grepl("BLACK|CARIBBEAN", ethnicity) ~ "black",
          grepl("HISPANIC|SOUTH AMERICAN", ethnicity) ~ "hispanic",
          grepl("HAWAIIAN", ethnicity) ~ "hawaiian",
          grepl("WHITE|PORTUGUESE", ethnicity) ~ "white",
          grepl(
            "OTHER|UNKNOWN|UNABLE|MULTI|DECLINED",
            ethnicity
          ) ~ "unknown"
        ),
      ethnicity = NULL,
      # Get 30-day mortality
      mortality_30d = ifelse(as.Date(dod) - as.Date(intime) <= 30 &
        !is.na(dod),
      TRUE, FALSE
      ),
      # Dichotomize fluid intake to > or <= 6L
      # in the first 24 hours of admission
      fluids_lt6L_24h = fluids_24h_l <= 6,
      fluids_24h_l = NULL,
      # Encode diabetes as factor
      diabetes = as.factor(diabetes),
      # Encode first_service as factor
      first_service = as.factor(first_service),
      # Encode qsofa and sirs as factors
      qsofa = as.factor(qsofa),
      sirs = as.factor(sirs)
    ) %>%
    arrange(intime) %>%
    as.data.frame()
}


sepsis <- fread("./data/sepsis_patients.csv",
  na.strings = "NULL"
) %>%
  clean_sepsis_data()

y <- as.numeric(sepsis$mortality_30d)
X <- model.matrix(
  ~ age + gender + diabetes +
    elixhauser_hospital + sofa +
    qsofa,
  data = sepsis
) %>%
  as.data.frame() %>%
  mutate(`(Intercept)` = NULL)

treatment <- as.numeric(sepsis$fluids_lt6L_24h)
SL.library <- c(
  "SL.earth",
  "SL.gam",
  "SL.glmnet",
  "SL.glm.interaction",
  "SL.ranger",
  "SL.xgboost"
)

# times at which to compute the confidence sequences
times <- unique(round(logseq(1000, nrow(sepsis), n = 30)))
t_opt <- 1500
alpha <- 0.1


cs_unadj <-
  asymptotic_confseq(
    x = (treatment == 1) * y / cumul_mean(treatment, regularizer_obs = 1) -
      (treatment == 0) * y / (cumul_mean(1 - treatment, regularizer_obs = 1)),
    t_opt = t_opt,
    alpha = alpha
  )

cs_unadj$l <- cs_unadj$l[times]
cs_unadj$u <- cs_unadj$u[times]

cs_glm <- confseq_ate(
  y = y,
  X = X,
  treatment = treatment,
  regression_fn_1 =
    get_SL_fn(
      SL.library = c("SL.glm"),
      family = binomial()
    ),
  propensity_score_fn =
    get_SL_fn(
      SL.library = c("SL.glm"),
      family = binomial()
    ),
  t_opt = t_opt,
  alpha = alpha,
  times = times,
  n_cores = 8,
  cross_fit = TRUE
)

cs_SL <- confseq_ate(
  y = y,
  X = X,
  treatment = treatment,
  regression_fn_1 = get_SL_fn(
    SL.library = SL.library,
    family = binomial()
  ),
  propensity_score_fn = get_SL_fn(
    SL.library = SL.library,
    family = binomial()
  ),
  t_opt = t_opt,
  alpha = alpha,
  times = times,
  n_cores = 8,
  cross_fit = TRUE
)

r_data_dir <- "./"
save(cs_unadj,
  cs_glm,
  cs_SL,
  times,
  file = paste(r_data_dir, "sepsis_cs.RData", sep = "")
)
