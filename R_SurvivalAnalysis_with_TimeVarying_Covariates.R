###############################################################################
# Created by: Ali Mirzazadeh
# For questions or assistance: ali.mirzazadeh@ucsf.edu
#
# Addicts Survival Analysis & Heart Transplant Survival Analysis
###############################################################################

# -------------------------
# 0. Setup
# -------------------------
# Install packages if needed:
# install.packages(c("tidyverse","readxl","survival","survminer","broom"))

library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(broom)

# ---- Change this to the folder that holds your data ----
setwd("/Users/alimirzazadeh1/Documents/GitHub/survival")

# -------------------------------------------------------

# ------------------------------------------------------------------------------
# ======================== Addicts Survival Analysis ===========================
# ------------------------------------------------------------------------------

# 1. Load data ---------------------------------------------------------------
addicts <- read_excel("ADDICTS.xlsx", sheet = "data") %>% as_tibble()

# Quick look
glimpse(addicts)
summary(addicts)
# Show first 10 rows
print(addicts %>% slice_head(n = 10))

# Note / assumptions about variables (adapt if your variable names differ):
# survt  = follow-up time
# status = event indicator (1 = relapse/event, 0 = censored)
# id     = subject identifier
# clinic = clinic indicator (1/2 in input) -> convert to 0/1
# prison = history of imprisonment (assumed 0/1)
# dose   = methadone dose (numeric)

# Recode clinic so it is 0/1 if it was 1/2
if("clinic" %in% names(addicts)) {
  addicts <- addicts %>% mutate(clinic = clinic - 1)
}

# 2. Declare survival object -------------------------------------------------
# For many analyses we'll use the counting-process style Surv(start, stop, event)
# For simple Kaplan-Meier and basic Cox we can use Surv(survt, status)

addicts <- addicts %>%
  mutate(
    survt = as.numeric(survt),
    status = as.integer(status)
  )

surv_obj <- Surv(time = addicts$survt, event = addicts$status)

# 3. Kaplan–Meier estimates and log-rank test --------------------------------

# Overall KM
fit_km <- survfit(surv_obj ~ 1, data = addicts)
print(fit_km)
ggsurvplot(fit_km, conf.int = TRUE,
           title = "KM: All subjects",
           xlab = "Time",
           ylab = "Survival probability")

# KM by clinic
fit_km_clinic <- survfit(surv_obj ~ clinic, data = addicts)
print(fit_km_clinic)
ggsurvplot(fit_km_clinic, conf.int = TRUE,
           title = "KM by clinic",
           xlab = "Time",
           ylab = "Survival probability",
           legend.title = "clinic")

# Log-rank test
logrank <- survdiff(surv_obj ~ clinic, data = addicts)
logrank
# print p-value
1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)

# 4. Cox proportional hazards model (basic PH) -------------------------------
cox_basic <- coxph(Surv(survt, status) ~ prison + clinic + dose, data = addicts)
summary(cox_basic)

# 5. Test proportional hazards assumption ------------------------------------

ph_test <- cox.zph(cox_basic)
ph_test
plot(ph_test) # visually inspect; lines = scaled schoenfeld residuals vs time


# 6. Heaviside (step-function) time-dependent effect -------------------------
# Create time-dependent covariate that equals clinic after 365 days (counting-process data)
# Approach: convert to counting-process dataset using survSplit
# We'll split at day 365 (and at event times if desired)

addicts_split <- survSplit(data = addicts,
                           cut = c(365), # cutoff at 365 days
                           end = "survt",
                           event = "status",
                           start = "tstart",
                           id = "id_split") %>%
  mutate(clinic_after365 = clinic * (tstart >= 365 | (tstart < 365 & survt >= 365)))

# Fit Cox with clinic baseline + clinic after 365 (Heaviside)
cox_heaviside <- coxph(Surv(tstart, survt, status) ~ prison + dose + clinic + clinic_after365, data = addicts_split)
summary(cox_heaviside)

# Alternatively use tt() with step function:
m2_step <- coxph(Surv(survt, status) ~ prison + dose + clinic + tt(clinic),
                 data = addicts,
                 tt = function(x, t, ...){ x * as.numeric(t >= 365) })
summary(m2_step)

# 7. Manual vs automatic time-varying coefficient (linear time interaction)
# Manual approach (must use counting-process data to have time in each row)
# Use survSplit to cut the dataset finely (e.g., at each event time) or keep as is but use tstop for time
# Here we'll use tstop = survt and tstart = 0 for each subject (no split)
addicts_counting <- addicts %>% mutate(tstart = 0, tstop = survt)

# Create clinic * time (using tstop to represent time at end of interval)
addicts_counting <- addicts_counting %>%
  mutate(clinic_t = clinic * tstop)

# Cox with explicit interaction term (interprets clinic_t as linear interaction)
cox_manual_interact <- coxph(Surv(tstart, tstop, status) ~ prison + dose + clinic + clinic_t, data = addicts_counting)
summary(cox_manual_interact)

# Equivalent using tt() (automatic)
cox_ttv_manual <- coxph(Surv(survt, status) ~ prison + dose + clinic + tt(clinic),
                        data = addicts,
                        tt = function(x, t, ...){ x * t })
summary(cox_ttv_manual)

# NOTE:
# - The manual approach multiplies clinic by the chosen time (tstop here).
# - To be equivalent to tvc with tt(...=t), you must ensure the same time variable is used
#   and the counting-process structure is consistent.

# ------------------------------------------------------------------------------
# ===================== Heart Transplant Survival Analysis ======================
# ------------------------------------------------------------------------------

# -------------------------
# 0. Setup
# -------------------------
# Install packages if needed:
# install.packages(c("tidyverse","readxl","survival","survminer","broom"))

library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(broom)

# ---- Change this to the folder that holds your data ----
setwd("/Users/alimirzazadeh1/Documents/GitHub/survival")


# 1. Load data ----------------------------------------------------------------
heart <- read_excel("HEARTT.xlsx", sheet = "data") %>% as_tibble()
glimpse(heart)
summary(heart)
print(heart %>% slice_head(n = 10))

# Assumed variables:
# stime = follow-up time
# died  = event indicator (1=died, 0=censored)
# id    = subject id
# wait  = time to transplant (0 or NA for never transplanted)
# age, surgery, year = covariates

# 2. Naive (baseline) Cox model — treating transplant as fixed -----------------
heart <- heart %>%
  mutate(stime = as.numeric(stime),
         died = as.integer(died))

# KM by transplant status (treated as baseline)
fit_km_ht <- survfit(Surv(stime, died) ~ transplant, data = heart)
ggsurvplot(fit_km_ht, conf.int = TRUE,
           title = "Heart transplant: KM by transplant (naive)")

# Naive Cox treating transplant as baseline covariate (this is biased if transplant occurs during follow-up)
cox_naive_transplant <- coxph(Surv(stime, died) ~ transplant, data = heart)
summary(cox_naive_transplant)

# 3. Time-dependent transplant variable (correct approach) ---------------------
# Use tmerge to create time-dependent transplant indicator correctly.
# Prepare baseline data for tmerge. tmerge will create counting-process style data.

# create a baseline dataset with start = 0, stop = stime, event = died
heart_base <- heart %>%
  mutate(tstart = 0, tstop = stime, event = died)

# If wait==0 or missing represent no transplant, set NA for transplant time; else use wait
# tmerge needs the time of the (time-dependent) transplant event.
heart_base <- heart_base %>%
  mutate(trans_time = ifelse(is.na(wait) | wait == 0, NA, wait))

# tmerge: add a time-dependent variable 'posttran' that becomes 1 at trans_time
# We'll use tmerge from survival package:
heart_tm <- tmerge(
  data1 = heart_base,
  data2 = heart_base,
  id = id,
  tstart = tstart,
  tstop = tstop,
  event = event(stime, event) # keep original event
)

# Add transplant as a time-dependent single event (0->1 at trans_time)
# tdc() uses the 'trans_time' column to add an indicator that switches at that time
heart_tm <- tmerge(heart_tm, heart_base, id = id, posttran = tdc(trans_time))

# Inspect resulting counting-process data
heart_tm %>% group_by(id) %>% summarize(nrows = n(), .groups = "drop") %>% print(n = Inf)

# Fit Cox with time-dependent posttran (counting-process style)
cox_td_transplant <- coxph(Surv(tstart, tstop, event) ~ age + posttran + surgery + year, data = heart_tm)
summary(cox_td_transplant)

# This is equivalent to Stata stsplit + posttran generation approach:
# - Before transplant: posttran == 0
# - After transplant:  posttran == 1
# - Never transplanted: posttran stays 0 for all intervals

