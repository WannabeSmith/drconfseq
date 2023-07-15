library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

clean_sepsis_data <- function(data) {
  data %>%
    lazy_dt(immutable = TRUE) %>%
    mutate(
      # Simplifying ethnicity column according to NIH guidelines
      # https://grants.nih.gov/grants/guide/notice-files/not-od-15-089.html
      ethnicity_NIH =
        case_when(
          grepl('AMERICAN INDIAN', ethnicity) ~ 'native',
          grepl('ASIAN|MIDDLE EASTERN', ethnicity) ~ 'asian',
          grepl('BLACK|CARIBBEAN', ethnicity) ~ 'black',
          grepl('HISPANIC|SOUTH AMERICAN', ethnicity) ~ 'hispanic',
          grepl('HAWAIIAN', ethnicity) ~ 'hawaiian',
          grepl('WHITE|PORTUGUESE', ethnicity) ~ 'white',
          grepl('OTHER|UNKNOWN|UNABLE|MULTI|DECLINED',
                ethnicity) ~ 'unknown'
        ),
      ethnicity = NULL,
      # Get 30-day mortality
      mortality_30d = ifelse(as.Date(dod) - as.Date(intime) <= 30 &
                               !is.na(dod),
                             TRUE, FALSE),
      # Dichotomize fluid intake to > or <= 6L
      # in the first 24 hours of admission
      fluids_lt6L_24h = fluids_24h_l <= 6,
      fluids_24h_l,
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
