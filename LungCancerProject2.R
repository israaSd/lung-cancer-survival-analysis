library(readxl)
library(MASS)
library(survival)
library(survminer)

# read table

data = read_excel("DataProject.xlsx")
dt = read.csv("doublingTime.csv")
dt$DT = dt$DT* 365.25
summary(dt)
hist(dt$DT,  breaks = 60,col="light blue", main="doubling time", xlab="doubling time")


#Data Cleaning:
####################################################################################
# Remove missing values
data_clean <- data[!is.na(data$TumorSize), ]
summary((data_clean))

set.seed(123)
M = nrow(data)
# fitting doubling time 
# gamma
# fitdistr: Maximum-likelihood fitting of univariate distributions, allowing parameters to be held fixed if desired.
X1 = fitdistr(dt$DT, "gamma")
X1

#Kolmogorov-Smirnov (K-S) test [Dn =supx ∣Fn(x)−F(x)∣] is a non-parametric test employed to check whether the probability distributions of a sample and a control distribution, or two samples are equal
ks.test(dt$DT,"pgamma", 1.962174675 ,  0.015487722)

#logLik computes the log-likelihood of observing the training data given the model -- this is precisely the parameter that was maximised during model fitting.
logLik(X1)  #-745.0714

#WEIBULL
X2<-fitdistr(dt$DT , "weibull")
X2
ks.test(dt$DT,"pweibull", 1.43849737 ,  140.35103769)
logLik(X2) #-747.0852

#Exp
X3<-fitdistr(dt$DT , "exponential")
X3
ks.test(dt$DT,"pexp", 0.0078931390)
logLik(X3) #-759.429

#log normal
X4<-fitdistr(dt$DT, "lognormal")
X4
ks.test(dt$DT,"plnorm", 4.56580355  ,  0.77142052 )
logLik(X4) #-744.2787

set.seed(123)
DT <- rlnorm(M,  4.56580355  , 0.77142052 )
for (i in 1:M) {
  if (DT [i] > 440) {DT[i]=mean(DT)}
  else {DT[i]=DT[i]}}
summary(DT) 


# Lead time = how much earlier screening finds the cancer compared with clinical detection.
# sojourn time ~ exponential(λ)

# Screened patients
screened = subset(data_clean, Screening == "Screened")
mean(screened$TumorSize)
s1 = fitdistr(screened$TumorSize, "gamma")
s1

ks.test(screened$TumorSize,"pgamma", 2.902469419 ,  0.074598004 )
logLik(s1) # -818.126

s2 = fitdistr(screened$TumorSize, "weibull")
s2
ks.test(screened$TumorSize,"pweibull", 1.51603448 ,  43.57489422  )
logLik(s2) #-834.5053 

s3 = fitdistr(screened$TumorSize, "exponential")
s3
ks.test(screened$TumorSize,"pexp", 0.025701584  )
logLik(s3) #-862.3225

s4 = fitdistr(screened$TumorSize, "lognormal")
s4
ks.test(screened$TumorSize,"plnorm", 3.47915447 ,  0.58899900   )
logLik(s4) # -808.221

set.seed(123)
D0<-rlnorm(M, 3.47915447 ,  0.58899900)
for (i in 1:M) {
  if (D0[i]> 100||D0[i]< 10) {D0[i]=mean(D0)}
  else {D0[i]=D0[i]}}
summary(D0)

# unscreened patients

unscreened = subset(data_clean, Screening == "Unscreened")
mean(unscreened$TumorSize)

us1 = fitdistr(unscreened$TumorSize, "gamma")
us1

ks.test(unscreened$TumorSize,"pgamma",2.122919387 ,  0.032222334  )
logLik(us1) # -2977.072

us2 = fitdistr(unscreened$TumorSize, "weibull")
us2
ks.test(unscreened$TumorSize,"pweibull", 1.34934787  , 72.55716866 )
logLik(us2) # -3008.081

us3 = fitdistr(unscreened$TumorSize, "exponential")
us3
ks.test(unscreened$TumorSize,"pexp", 0.0151772830 )
logLik(us3) # -3055.706 

us4 = fitdistr(unscreened$TumorSize, "lognormal")
us4
ks.test(unscreened$TumorSize,"plnorm", 3.93417587 ,  0.70981514  )
logLik(us4) # -2951.104

set.seed(123)
D1<-rlnorm(M,  3.93417587 ,  0.70981514 )
for (i in 1:M) {
  if (D1[i]> 100||D1[i]< 10) {D1[i]=mean(D0)}
  else {D1[i]=D1[i]}}
summary(D1)


### Sojourn Time #####
set.seed(123)
library(triangle)
dt0 = rtriangle(M, a = 1, b = 1.4, c = 1.2)

set.seed(123)
SojournTime = 3 * DT * (log(D1 / dt0) / log(2))
summary(SojournTime)
SojournTime = replace (SojournTime, SojournTime<0,mean(SojournTime))
SojournTime = replace (SojournTime, SojournTime>3650,mean(SojournTime))
summary(SojournTime)

hist(SojournTime, pch=16,col="lightblue", xlab="Sojourn Time (Days)")

MST = mean(SojournTime,na.rm=T) # 1626.527

fitting1 = fitdistr(SojournTime,"lognormal")
fitting1
ks.test(SojournTime,"plnorm",7.17723201 ,  0.72709325)
logLik(fitting1) # -7400.058

fitting2 = fitdistr(SojournTime,"gamma")
fitting2
ks.test(SojournTime,"pgamma", 2.458312e+00 ,  1.511388e-03  )
logLik(fitting2) #  -7342.899 

fitting3 = fitdistr(SojournTime,"exponential")
fitting3
ks.test(SojournTime,"pexp", 6.148069e-04 )
logLik(fitting3) #-7504.417

fitting4 = fitdistr(SojournTime,"weibull")
fitting4
ks.test(SojournTime,"pweibull",  )
logLik(fitting4)

# Transition rate λ = 1 / MST
lambda = 1 / MST # 0.0006148069

##################### Lead Time

#############################LeadTime Schwartz##################################
set.seed(123)
leadtime = 3*DT*((log(D1/D0))/log(2))
leadtime = replace (leadtime, leadtime<0,mean(leadtime))
leadtime = replace (leadtime,leadtime> SojournTime,mean(leadtime))
summary(leadtime, na.rm=T)
hist(leadtime, pch=16,col="lightblue", xlab="lead timeTime")

t.test(leadtime, na.rm=T)$conf.int

quantile(leadtime, na.rm=T)

CorrectedSurvival = data$SurvivalInDays-leadtime
mean(CorrectedSurvival,na.rm=T)
t.test(CorrectedSurvival,na.rm=T)$conf.int
mean(data$SurvivalInDays,na.rm=T) - mean(CorrectedSurvival,na.rm=T)

data$CorrectedSurvival = data$SurvivalInDays-leadtime
summary(data)
data$CorrectedSurvival[data$CorrectedSurvival < 0] = median(data$CorrectedSurvival)
summary(data$CorrectedSurvival)
###################################### kaplan Mair & log-rank test

# Recode DeathStatus into binary
data$DeathStatus_num <- ifelse(data$DeathStatus == "Death", 1, 0)
# data$CorrectedSurvival

# Factor : Gender
# Kaplan-Meier fit
fit_G = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ Gender, data = data)
fit_G1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ Gender, data = data)
# Plot KM curves
ggsurvplot(fit_G, data = data, pval = TRUE) # curves overlap : likely no effect
ggsurvplot(fit_G1, data = data, pval = TRUE) # # Male lower than Female

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ Gender, data = data) 
# p = 0.05 no evidence the factor gender affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ Gender, data = data) 
# p = 0.01 significantly affects 

# Factor : AgeGroup
# Kaplan-Meier fit
fit_AG = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ AgeGroup, data = data)
fit_AG1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ AgeGroup, data = data)
# Plot KM curves
ggsurvplot(fit_AG, data = data, pval = TRUE) # curves overlap : likely no effect
ggsurvplot(fit_AG1, data = data, pval = TRUE) # curves overlap : likely no effect

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ AgeGroup, data = data) 
# p = 0.05 no evidence the factor affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ AgeGroup, data = data)
# p = 0.7 no evidence the factor affects survival

# Factor : Smoker
# Kaplan-Meier fit
fit_S = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ Smoker, data = data)
fit_S1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ Smoker, data = data)

# Plot KM curves
ggsurvplot(fit_S, data = data, pval = TRUE) # smoker (reach events (death) earlier) lower than No smoker (survive longer) : effect
# curves look separated : likely the factor affects survival
ggsurvplot(fit_S1, data = data, pval = TRUE) # curves look separated : likely the factor affects survival

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ Smoker, data = data)
# p = p= 7e-05  factor smoker significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ Smoker, data = data)
# p = p= 0.004  factor smoker significantly affects survival

# AbstinenceStatus
fit_Abst = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ AbstinenceStatus, data = data)
fit_Abst1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ AbstinenceStatus, data = data)

# Plot KM curves
ggsurvplot(fit_Abst, data = data, pval = TRUE) # Non Abstinet Smokerrelate lower than Abstinet Smokerrelated : effect
# curves look separated : likely the factor affects survival
ggsurvplot(fit_Abst1, data = data, pval = TRUE)
# curves look separated : likely the factor affects survival

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ AbstinenceStatus, data = data)
# p = 0.01 significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ AbstinenceStatus, data = data)
# p = 0.002 significantly affects survival

# HBV
fit_HBV = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ HBV, data = data)
fit_HBV1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ HBV, data = data)

# Plot KM curves
ggsurvplot(fit_HBV, data = data, pval = TRUE) # curves overlap : likely no effect
ggsurvplot(fit_HBV1, data = data, pval = TRUE) # curves overlap : likely no effect

# Log-rank test
survdiff(Surv(SurvivalInDays, data$DeathStatus_num) ~ HBV, data = data)
# p=0.8 no evidence the factor affects survival
survdiff(Surv(CorrectedSurvival, data$DeathStatus_num) ~ HBV, data = data)
# p=0.8 no evidence the factor affects survival

# HIV 
fit_HIV = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ HIV, data = data)
fit_HIV1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ HIV, data = data)

# Plot KM curves
ggsurvplot(fit_HIV, data = data, pval = TRUE) # No HIV lower than HIV
ggsurvplot(fit_HIV1, data = data, pval = TRUE) # No HIV lower than HIV

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ HIV, data = data)
# p = 0.002 significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ HIV, data = data)
# p = 0.05

# Diabetes
fit_d = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ Diabetes, data = data)
fit_d1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ Diabetes, data = data)

# Plot KM curves
ggsurvplot(fit_d, data = data, pval = TRUE) # curves overlap : likely no effect
ggsurvplot(fit_d1, data = data, pval = TRUE) # a little more separate than fit_d

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ Diabetes, data = data)
# p= 0.3 no evidence the factor affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ Diabetes, data = data)
# p = 0.09 

# Alcohol_consumption
fit_ALC = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ Alcohol_consumption, data = data)
fit_ALC1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ Alcohol_consumption, data = data)

# Plot KM curves
ggsurvplot(fit_ALC, data = data, pval = TRUE) # curves overlap : likely no effect
ggsurvplot(fit_ALC1, data = data, pval = TRUE) # a little more separate

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ Alcohol_consumption, data = data)
# p = 0.6 no evidence the factor affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ Alcohol_consumption, data = data)
# p = 0.3

# Thrombosis 
fit_Thr = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ Thrombosis, data = data)
fit_Thr1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ Thrombosis, data = data)

# Plot KM curves
ggsurvplot(fit_Thr, data = data, pval = TRUE) # Thrombosis lower than No thromb. : effect
ggsurvplot(fit_Thr1, data = data, pval = TRUE) # Thrombosis lower than No thromb. : effect

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ Thrombosis, data = data)
# p= <2e-16 significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ Thrombosis, data = data)
# p= 5e-13 significantly affects survival

# DifuseCancer
fit_dif = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ DifuseCancer, data = data)
fit_dif1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ DifuseCancer, data = data)

# Plot KM curves
ggsurvplot(fit_dif, data = data, pval = TRUE) # DifuseCancer lower than No : effect
ggsurvplot(fit_dif1, data = data, pval = TRUE) # DifuseCancer lower than No : effect

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ DifuseCancer, data = data)
# p= <2e-16  significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ DifuseCancer, data = data)
# p= <1e-16  significantly affects survival

# MetastaticCancer 
fit_metast = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ MetastaticCancer, data = data)
fit_metast1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ MetastaticCancer, data = data)

# Plot KM curves
ggsurvplot(fit_metast, data = data, pval = TRUE) # MetastaticCancer lower than No : effect
ggsurvplot(fit_metast1, data = data, pval = TRUE) # MetastaticCancer lower than No : effect

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ MetastaticCancer, data = data)
# # p= <2e-16  significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ MetastaticCancer, data = data)
# p= <2e-10  significantly affects survival

# CurativeTreatment
fit_cur = survfit(Surv(SurvivalInDays, DeathStatus_num) ~ CurativeTreatment, data = data)
fit_cur1 = survfit(Surv(CorrectedSurvival, DeathStatus_num) ~ CurativeTreatment, data = data)

# Plot KM curves
ggsurvplot(fit_cur, data = data, pval = TRUE) # Treatment without a curat lower than with : effect
ggsurvplot(fit_cur1, data = data, pval = TRUE) # Treatment without a curat lower than with : effect

# Log-rank test
survdiff(Surv(SurvivalInDays, DeathStatus_num) ~ CurativeTreatment, data = data)
# p= <2e-16  significantly affects survival
survdiff(Surv(CorrectedSurvival, DeathStatus_num) ~ CurativeTreatment, data = data)
# p= <2e-16  significantly affects survival

############# COX Regression

cox_obs = coxph(Surv(SurvivalInDays, DeathStatus_num) ~ Gender + Age + AgeGroup + Smoker + HBV + HIV + TumorSize + Diabetes + Alcohol_consumption + Thrombosis + DifuseCancer + MetastaticCancer + CurativeTreatment +Treatment1 +Treatment2 +Treatment3 +Treatment4 +Treatment5 +Treatment6, data = data)

cox_corr = coxph(Surv(CorrectedSurvival, DeathStatus_num) ~ Gender + Age + AgeGroup + Smoker + HBV + HIV + TumorSize + Diabetes + Alcohol_consumption + Thrombosis + DifuseCancer + MetastaticCancer + CurativeTreatment +Treatment1 +Treatment2 +Treatment3 +Treatment4 +Treatment5 +Treatment6, data = data)

# exp(coef) : Hazard Ratio (HR)
# p : significance

capture.output(coxph(Surv(SurvivalInDays, data$DeathStatus_num) ~ Gender + Age + AgeGroup + Smoker + HBV + HIV + TumorSize + Diabetes + Alcohol_consumption + Thrombosis + DifuseCancer + MetastaticCancer + CurativeTreatment +Treatment1 +Treatment2 +Treatment3 +Treatment4 +Treatment5 +Treatment6, , data = data), coxph(Surv(CorrectedSurvival, data$DeathStatus_num) ~ Gender + Age + AgeGroup + Smoker + HBV + HIV + TumorSize + Diabetes + Alcohol_consumption + Thrombosis + DifuseCancer + MetastaticCancer + CurativeTreatment +Treatment1 +Treatment2 +Treatment3 +Treatment4 +Treatment5 +Treatment6, data = data), file = "Cox_model.txt")
################################################3

library(broom)
library(dplyr)

cox_table_obs <- tidy(cox_obs, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    HR = estimate,
    CI_lower = conf.low,
    CI_upper = conf.high,
    p_value = p.value
  )

cox_table_corr <- tidy(cox_corr, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  rename(
    HR = estimate,
    CI_lower = conf.low,
    CI_upper = conf.high,
    p_value = p.value
  )

library(knitr)
library(kableExtra)

kable(
  cox_table_obs,
  caption = "Appendix Table A1. Multivariable Cox regression using observed survival time",
  digits = 3
) %>%
  kable_styling(full_width = FALSE)

kable(
  cox_table_corr,
  caption = "Appendix Table A2. Multivariable Cox regression using lead-time–corrected survival",
  digits = 3
) %>%
  kable_styling(full_width = FALSE)
