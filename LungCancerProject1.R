setwd('C:\\Users\\admin\\Desktop\\stat 504\\project 1')

# read table

library(readxl)
data = read_excel("DataProject.xlsx")

# Select categorical and continuous variables :

#categorical_vars <- c("gender","age_group","smoker","screened","diabetes","hbv","hcv","nash","alcohol","metastases","curative","death")
#continuous_vars <- c("age","tumor_size","surv_days")

# summary data

library(summarytools)
dfSummary(data)
 # Gender 894: female 156 (17.4%)  ,  MAle 738 (82.6%)
 # Age 891 (99.7%): mean(sd): 67.4 (11.3) , min: 17 , med: 68, max: 99
 # Age Group 891: <= 70 : 504 (56.6%) , > 70 : 387 (43.4%)
 # Smoker 894: No Smoker related Cancer : 312 (34.9%) , Smoker related cancer :  582 (65.1%)
 # AbstinenceStatus 549 (61.4%) : Abstinet Smokerrelated Ca : 305 (55.6%) , Non Abstinet Smokerrelate : 244 (44.4%)
 # HBV 795 (88.9%):  No : 723 (90.9%) , Yes : 72 (9.1%)
 # HIV 800 (89.5%): No :  698 (87.2%) , Yes :  102 (12.8%)
 # other ?
 # Screening 884 (98.9%) : Screened : 199 (22.5%) , Unscreened  : 685 (77.5%)
 # TumorSize in mm 113 (87.4%): mean (sd) : 59.4 (56.4) , min : 2, med: 45, max: 900
 # Diabetes 885 (99%): NO 595 (66.9%) , Yes 293 (33.1%)
 # Alcohol_consumption 852 (95.3%): No 551 (64.7%), Yes 251 (35.3%)
 # Thrombosis 881 (98.5%) : No 630 (71.5%) , Yes 251 (28.5%)
 # Criteria 891 (99.7%) : In (Tumor Size<5cm) 254 (28.5%) , Out (>5cm) 637 (71.5%)
 # Cancer Stages 740 (82.8%) : A (early stage) 113 (15.3%), B (intermediate stage) 41 (5.5%)
 #                             C (Advance stage) 346 (46.8%), D (terminal stage) 240 (32.4%)
 # DifuseCancer 885 (99%): No 745 (84.2%), Yes 113 (12.7%)
 # MetastaticCancer 887 (99.2%): No 774 (87.3%), Yes 113 (12.7%)
 # CurativeTreatment 865 (96.8%): Treatment with a curative 173 (20.0%), Treatment without a curat 692 (80.0%)
 # Treatment1 875 (97.9%): No 787 (89.9%), Yes 88 (10.1%)
 # Treatment2 871 (97.4%): No 800 (91.8%), Yes 71 ( 8.2%)
 # Treatment3 875 (97.9%): No 729 (83.3%), Yes 146 (16.7%)
 # Treatment4 871 (97.4%): No 718 (82.4%), Yes 153 (17.6%)
 # Treatment5 871 (97.4%): No 542 (62.2%), Yes 329 (37.8%)
 # Treatment6 872 (97.5%): No 803 (92.1%), Yes 69 ( 7.9%)
 # DeathStatus 818 (91.5%): Death 601 (73.5%), No 217 (26.5%)
 # SurvivalInDays 108 (87.9%): mean (sd) : 328.1 (389.1), min:1, med: 209.5, max:2155

# save summary as text file :
capture.output(dfSummary(data), file = "summary_data.txt")


# Histogram of tumor size:
library(ggplot2)

ggplot(data, aes(x = TumorSize)) +
  geom_histogram(binwidth = 10, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Tumor Size", x = "Tumor size (mm)", y = "Count")

# Boxplot of survival
ggplot(data, aes(y = SurvivalInDays)) +
       geom_boxplot(fill = "lightgreen") +
       labs(title = "Boxplot of Survival (days)", y = "Survival time (days)")

# Mortality rate (and 95% CI)
library(epitools)
deaths = sum(data$DeathStatus == "Death", na.rm = TRUE)
n=sum(!is.na(data$DeathStatus))
prop.test(deaths, n) # rate= 0.7347188 , 95% CI: [0.7027879 ; 0.7644113]

# Rate of screened cancers
m=sum(!is.na(data$Screening))
screened = sum(data$Screening == "Screened", na.rm = TRUE)
prop.test(screened, m)  # rate= 0.2251131, 95% CI [0.1982654 ; 0.2543883]

# Risk of death and RR for smokers vs non-smokers
table(data$Smoker, data$DeathStatus)
# risk of death for smoker
RS  = 418/582
# Risk of death for no smoker
RNS = 183/312
# relative risk of death for smoker
RRS = RS/RNS # the risk of cancer death is increased by 22.449% for smokers compared with no smokers

# Risk of death and RR for Gender
table(data$Gender, data$DeathStatus)
# risk of death for Female
RF = 93/156
# risk of death for MAle
RM = 508/738
# relative risk of death for male
RRG = RM/RF # the risk of cancer death is 15.46% more for male than female

# Risk of death and RR for AgeGroup
table(data$AgeGroup, data$DeathStatus)
# risk of death for less or = 70 
R1 = 325/504
# risk of death for more than 70
R2 = 275/387
# relative risk of death for AgeGroup
RRAG = R2/R1   # the risk of cancer death is 10.197% more for people age>70 than age<=70

# Risk of death and RR for HBV
table(data$HBV, data$DeathStatus)
# risk of death for HBV
RHBV = 42/72
# risk of death for NO HBV
RNHBV = 484/723
# relative risk of death for HBV
RRHBV = RNHBV/RHBV # risk of cancer death is 14.7599% more for people with HBV than with HBV

# Risk of death and RR for HIV
table(data$HIV, data$DeathStatus)
# risk of death for HIV
RHIV = 55/102
# risk of death for NO HIV
RNHIV = 476/698
# relative risk of death for HIV
RRHIV = RRNHIV/RRHIV # risk of cancer death is 26.47% more for people with HIV than with HIV

# risk of death for Screening
table(data$Screening, data$DeathStatus)
# risk of death for screened people
RSC = 110/199
# risk of death for unscreened people
RUNSC = 484/685
# relative risk of death for screening
RRSC = RUNSC/RSC # risk of death 27.82% more for unscreened people

# risk of death for Diabetes
table(data$Diabetes, data$DeathStatus)
# risk of death for Diabetes
RD = 199/293
# risk of death for NO Diabetes
RND = 393/592
# relative risk of death for Diabetes
RRD = RD/RND # risk of death is 2.3% more for Diabetes people

# risk of death for Alcohol_consumption
table(data$Alcohol_consumption, data$DeathStatus)
# risk of death for Alcohol consumption
RAC = 206/301
# risk of death for no alcohol consumption
RNAC = 362/551
# relative risk of death for alcohol consumption
RRAC = RAC/RNAC # risk of death is 4.17% more for alcohol consumption

# risk of death for Thrombosis
table(data$Thrombosis, data$DeathStatus)
# risk of death for Thrombosis
RTH = 206/251
# risk of death for No Thrombosis
RNTH = 384/630
# relative risk of death for Thrombosis
RRTH = RTH/RNTH # 34.489% more for Thrombosis people

# risk of death for DifuseCancer
table(data$DifuseCancer, data$DeathStatus)
# risk of death for DifuseCancer
RDC = 119/140
# risk of death for No DifuseCancer
RNDC = 477/745
# relative risk of death for DifuseCancer
RRDC = RDC/RNDC # 32.75% more for DifuseCancer

# risk of death for MetastaticCancer
table(data$MetastaticCancer, data$DeathStatus)
# risk of death for MetastaticCancer
RMC = 91/113
# risk of death for No MetastaticCancer
RNMC = 504/774
# relative risk of death for MetatstaticCancer
RRMC = RMC/RNMC # 23.67% more for MetastaticCancer

# risk of death for CurativeTreatment 
table(data$CurativeTreatment, data$DeathStatus)
# risk of death for treat with curative intent
RCI = 70/173
# risk of death for treat without curative intent
RNCI = 518/692
# relative risk of death for treatment with curative intent
RRCI = RCI/RNCI # 85% more for treatment with curative intent

# impact of tumor size
# Recode DeathStatus into binary
data$DeathStatus_num <- ifelse(data$DeathStatus == "Death", 1, 0)
# Fit the logistic regression model
model = glm(DeathStatus_num ~ TumorSize, data = data, family = binomial)
summary(model)
exp(coef(model))  #1.458071    1.009914 
#For every 1 mm increase in tumor size, the odds of death increase by about 0.99%

#Fit a simple linear regression
# Example
model_linear = lm(SurvivalInDays ~ TumorSize, data = data)
summary(model_linear)
# Intercept (425.73) → predicted survival time (in days) when tumor size = 0 mm.
# Slope (−1.1643) → average change in survival days for each additional millimeter of tumor size.
# For every 1 mm increase in tumor size, the average survival time decreases by about 1.16 days.


# Survival analysis
# Lead Time Bias
hist(data$TumorSize,  breaks = 100,col="light blue", main="Tumor size of no screened patients", xlab="Tumor size (mm)")