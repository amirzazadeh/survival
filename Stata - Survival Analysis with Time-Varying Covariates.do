
* Created by Ali Mirzazadeh
* For questions or assistance, please contact:
* Email: ali.mirzazadeh@ucsf.edu

********************************************************************
* Addicts Survival Analysis
********************************************************************

clear all
set more off
* change the working directory to the folder that contains the data file on your computer
cd "/Users/alimirzazadeh1/Documents/GitHub/survival"

*---------------------------------------------------------------*
* 1. Load data
*---------------------------------------------------------------*
clear
import excel "ADDICTS.xlsx", sheet("data") firstrow

* Recode clinic so that it is 0/1 instead of 1/2
* This makes interpretation of hazard ratios easier
replace clinic = clinic - 1

describe
summarize


*---------------------------------------------------------------*
* 2. Declare survival-time data
*---------------------------------------------------------------*
* survt  = follow-up time
* status = event indicator (1 = relapse/event, 0 = censored)
* id     = subject identifier

stset survt, failure(status==1) id(id)

* Check survival-time settings
stset
stdescribe


*---------------------------------------------------------------*
* 3. Kaplan–Meier estimates and log-rank test
*---------------------------------------------------------------*

* Overall Kaplan–Meier survival estimates
sts list

* KM estimates by clinic group
sts list, by(clinic) compare at(0 20 to 1080)

* Plot KM survival curves by clinic
* Visual inspection of survival differences
sts graph, by(clinic)

* Log-rank test (nonparametric comparison of survival curves)
sts test clinic


*---------------------------------------------------------------*
* 4. Cox proportional hazards model (basic PH model)
*---------------------------------------------------------------*

* Fits standard Cox PH model:
* Assumes hazard ratios are constant over time
stcox prison clinic dose


*---------------------------------------------------------------*
* 5. Test proportional hazards assumption
*---------------------------------------------------------------*

* Fit standard Cox model
stcox prison clinic dose, cformat(%9.2f)

* Plot adjusted survival curves by clinic
* If curves cross, PH assumption may be violated
stcurve, survival at(clinic=0) at(clinic=1)

* Schoenfeld residual test of proportional hazards
* Tests whether effect changes over time
estat phtest, detail


* Likelihood ratio test for PH violation

* Model 1: Proportional hazards (clinic effect constant)
stcox prison clinic dose, cformat(%9.2f)
estimates store M1

* Model 2: Allow clinic effect to vary linearly with time
* tvc(clinic) creates clinic × time interaction
* texp(_t) specifies interaction with analysis time
stcox prison clinic dose, tvc(clinic) texp(_t) nohr
estimates store M2

* Likelihood ratio test:
* H0: clinic effect is proportional (no time interaction)
lrtest M1 M2 


*---------------------------------------------------------------*
* 6. Heaviside (step-function) time-dependent effect
*    Effect of clinic after 365 days
*---------------------------------------------------------------*

* This allows clinic to have one effect before 365 days
* and a different effect after 365 days
* (_t >= 365) is a step function

stcox prison dose clinic, cformat(%9.2f) tvc(clinic) texp(_t >= 365) nohr
stcox prison dose clinic, cformat(%9.2f) tvc(clinic) texp(_t >= 365)

* Interpretation:
* Baseline clinic coefficient = effect before 365 days
* Interaction coefficient     = additional effect after 365 days


*---------------------------------------------------------------*
* 7. Manual vs automatic time-varying coefficient
*    Linear time interaction for clinic
*---------------------------------------------------------------*

stset survt, failure(status==1) id(id)

* Manually create clinic × time interaction
cap drop clinic_t
generate clinic_t = clinic * _t

* Cox model including explicit interaction term
* Allows clinic hazard ratio to change linearly over time
stcox prison dose clinic clinic_t, cformat(%9.3f) nohr

* Equivalent model using tvc() option
* tvc(clinic) texp(_t) automatically creates clinic × time
stcox prison dose clinic, cformat(%9.3f) tvc(clinic) texp(_t) nohr

* If implemented correctly, both models should produce identical results.
* Differences usually indicate:
*   - clinic_t created before stset
*   - different time scale used
*   - data modified between models





********************************************************************
* Heart Transplant Survival Analysis
********************************************************************

clear all
set more off
* change the working directory to the folder that contains the data file on your computer
cd "/Users/alimirzazadeh1/Documents/GitHub/survival"


*---------------------------------------------------------------*
* 1. Load data
*---------------------------------------------------------------*
clear
import excel "HEARTT.xlsx", sheet("data") firstrow

describe
summarize
list in 1/10

*---------------------------------------------------------------*
* 2. Naïve (baseline) Cox model – treating transplant as fixed
*---------------------------------------------------------------*

* Declare survival-time data
stset stime, failure(died) id(id)

* Kaplan–Meier survival curves by transplant status
* (This treats transplant as if known at baseline)
sts graph, by(transplant)

* Cox model without hazard ratios (log-coefficients shown)
stcox transplant, cformat(%9.2f) nohr

* Cox model showing hazard ratios
stcox transplant, cformat(%9.2f)

* NOTE:
* This model treats transplant as a baseline covariate.
* It ignores the fact that transplant occurs during follow-up.
* This can introduce immortal time bias.


*---------------------------------------------------------------*
* 3. Time-dependent transplant variable (correct approach)
*---------------------------------------------------------------*

* Re-declare survival data (ensures clean setup)
stset stime, failure(died) id(id)

* Split follow-up time at each failure time
* This creates multiple time intervals per subject
stsplit, at(failures)

* Create a time-dependent transplant indicator:
* posttran = 1 only after transplant has occurred
* wait = time to transplant
* _t = current analysis time
cap drop posttran
generate posttran = wait < _t & wait != 0

* Recombine split data (clean episode structure)
stjoin

* Cox model with time-dependent transplant variable
stcox age posttran surgery year, cformat(%9.2f)

* NOTE:
* posttran now changes over time.
* Before transplant → posttran = 0
* After transplant  → posttran = 1
* Never transplanted → always 0
*
* This correctly models transplant as a time-varying exposure
* and removes immortal time bias.
