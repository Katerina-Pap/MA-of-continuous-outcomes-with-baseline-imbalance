#----------------------------------------------------------------------------------------------------------------------------
#     R code supplementing " The impact of trial baseline imbalances should be considered in systematic reviews: 
#     a methodological case study revisited
#     Author: Katerina Papadimitropoulou
#     Date:   November 2020
#----------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------
#                                   Important R libraries and 
#                              data as reported, incl. missing values 
#----------------------------------------------------------------------------------------------------------------------------
library(readxl)
library(reshape2)
library(lme4)
library(metafor)
library(nlme)

# Load datasets
data.AD <- read_excel("Trowman_withNAs.xlsx") # Load Trowman/calcium supplementation dataset
# data.AD <- read_excel("apnea_withNAs.xlsx") # Load apnea dataset

#--------------------------------------------------------------------------------------------------------------------------
#                       Start algebraic calculations and imputations of values  
#----------------------------------------------------------------------------------------------------------------------------

# Calculate post baseline mean from CFB and baseline 
data.AD$MeanFU  <- ifelse(is.na(data.AD$MeanFU), data.AD$MeanCFB + data.AD$MeanBaseline, data.AD$MeanFU)

# Calculate change score values from baseline and follow-up
data.AD$MeanCFB <-  ifelse(is.na(data.AD$MeanCFB), data.AD$MeanFU - data.AD$MeanBaseline, data.AD$MeanCFB)

## Calculate missing standard deviations from standard errors and vice versa

# Calcucate SD from SE 
data.AD$sdBaseline     <- ifelse(is.na(data.AD$sdBaseline), data.AD$seBaseline*sqrt(data.AD$NCFB), data.AD$sdBaseline)
data.AD$sdFU           <- ifelse(is.na(data.AD$sdFU), data.AD$seFU*sqrt(data.AD$NCFB), data.AD$sdFU)
data.AD$sdCFB          <- ifelse(is.na(data.AD$sdCFB), data.AD$seCFB*sqrt(data.AD$NCFB), data.AD$sdCFB)

# If not possible assume same SD at baseline and follow-up
data.AD$sdFU           <- ifelse(is.na(data.AD$sdFU), data.AD$sdBaseline, data.AD$sdFU)

# Calculate group correlations using Equation (B8) or impute from reported median correlations of the remaining studies
data.AD$Correlation    <- ifelse(is.na(data.AD$Correlation), (data.AD$sdBaseline^2+data.AD$sdFU^2-data.AD$sdCFB^2)/(2*data.AD$sdBaseline*data.AD$sdFU), 
                                    data.AD$Correlation)

data.AD$Correlation    <- ifelse(is.na(data.AD$Correlation) & (data.AD$group=="0"), tapply(data.AD$Correlation, data.AD$group, median, na.rm=T)[1],
                                    data.AD$Correlation)

data.AD$Correlation    <- ifelse(is.na(data.AD$Correlation) & (data.AD$group=="1"), tapply(data.AD$Correlation, data.AD$group, median, na.rm=T)[2], 
                                    data.AD$Correlation)
# Final calculations of SE from SD
data.AD$seBaseline     <- ifelse(is.na(data.AD$seBaseline), data.AD$sdBaseline/sqrt(data.AD$NCFB), data.AD$seBaseline)
data.AD$seFU           <- ifelse(is.na(data.AD$seFU), data.AD$sdFU/sqrt(data.AD$NCFB), data.AD$seFU)
data.AD$sdCFB          <- ifelse(is.na(data.AD$sdCFB), sqrt(data.AD$sdBaseline^2+data.AD$sdFU^2-2*data.AD$Correlation*data.AD$sdBaseline*data.AD$sdFU), data.AD$sdCFB)
data.AD$seCFB          <- ifelse(is.na(data.AD$seCFB), data.AD$sdCFB/sqrt(data.AD$NCFB), data.AD$seCFB)
data.AD
#----------------------------------------------------------------------------------------------------------------------------
#             Perform (standard) AD approaches, Follow-up analysis, Change scores analysis 
#                                  and recovering ANCOVA methods
#----------------------------------------------------------------------------------------------------------------------------

# For the first three methods the data need to be in wide format 
drop         <- which(colnames(data.AD) %in% "Study")
data.AD      <- data.AD[,-drop]
data.AD_wide <- dcast(melt(data.AD, id.vars=c("ID", "group")), ID~variable+group)

#----------------------------------------------------------------------------------------------
#                             Method 1:  Final scores
#----------------------------------------------------------------------------------------------

MA.random.final <- rma(m1i=MeanFU_1, m2i=MeanFU_0, sd1i=sdFU_1, sd2i=sdFU_0, n1i=NCFB_1, n2i=NCFB_0,
                       data=data.AD_wide, measure="MD", knha=TRUE)
summary(MA.random.final)

#----------------------------------------------------------------------------------------------
#                             Method 2:  Change scores
#----------------------------------------------------------------------------------------------

MA.random.changescores <- rma(m1i=MeanCFB_1, m2i=MeanCFB_0, sd1i=sdCFB_1, sd2i=sdCFB_0, n1i=NCFB_1, n2i=NCFB_0,
                              data=data.AD_wide, measure="MD", knha=TRUE)
summary(MA.random.changescores)

#----------------------------------------------------------------------------------------------
#                           Method 3:  Recovering ANCOVA estimates approach
#----------------------------------------------------------------------------------------------

# Calculate pooled standard deviations of baseline and follow-up values
sdpooledB <- with(data.AD_wide, sqrt((((NCFB_1 - 1)*(sdBaseline_1^2)) + (NCFB_0 - 1)*(sdBaseline_0^2))/((NCFB_1+NCFB_0)-2)))
sdpooledF <- with(data.AD_wide, sqrt((((NCFB_1 - 1)*(sdFU_1^2)) + (NCFB_0 - 1)*(sdFU_0^2))/((NCFB_1+NCFB_0)-2)))

# Calculate ancova estimate using formula from Senn et al. 2007
# using the pooled correlation 

ripooled       <- with(data.AD_wide, ((NCFB_1*Correlation_1*sdBaseline_1*sdFU_1 +  NCFB_0*Correlation_0 *sdBaseline_0*sdFU_0) ) 
                 /((NCFB_1+NCFB_0)*sdpooledB*sdpooledF))

ancova_est     <- with(data.AD_wide, (MeanFU_1-MeanFU_0)-ripooled*(sdpooledF/sdpooledB)*(MeanBaseline_1-MeanBaseline_0))

var_ancova_est <- with(data.AD_wide, sdpooledF^2*(1/NCFB_1)+sdpooledF^2*(1/NCFB_0))*(1-ripooled^2) # for different sample sizes from McKenzie and from Senn

se_ancovas_est <- with(data.AD_wide,sqrt(var_ancova_est))

MA.random.ANCOVA <- rma(yi=ancova_est , sei=se_ancovas_est , slab=data.AD$study, method="REML", knha=TRUE)
summary(MA.random.ANCOVA)

#----------------------------------------------------------------------------------------------
#                            Method 5:  Modified Trowman approach
#----------------------------------------------------------------------------------------------
diff            <- with(data.AD_wide, MeanBaseline_0-MeanBaseline_1)
modified.random <-  rma(m1i=MeanFU_1, m2i=MeanFU_0, sd1i=sdFU_1, sd2i=sdFU_0, n1i=NCFB_1, n2i=NCFB_0, measure="MD", 
                        mods=~diff, method="REML", data=data.AD_wide, knha=TRUE)
summary(modified.random)

# Modified Trowman with interaction

modified.random.INT <- rma(m1i=MeanFU_1, m2i=MeanFU_0, sd1i=sdFU_1, sd2i=sdFU_0, n1i=NCFB_1, n2i=NCFB_0, measure="MD", 
                           mods=~diff + MeanBaseline_1, method="REML", data=data.AD_wide, knha=TRUE)
summary(modified.random.INT)


# Same results as the main effect of modified Trowman approach can be obtained by using Fay-Herriot estimate of the small area estimation library
yi <- c(MA.random.final$yi[1:8])
vi <- MA.random.final$vi
require(sae) 
va <- eblupFH(formula = yi ~ diff, vardir = vi, method = "REML")

#----------------------------------------------------------------------------------------------------------------------------
#                                         Method 6:  Pseudo IPD approach
#                             First generate the pseudo IPD and then fit mixed effects models
#----------------------------------------------------------------------------------------------------------------------------

# Use data in long format: data.AD

# Generate the pseudo baselines and outcomes
data.IPD <- data.frame(study         = rep(data.AD$ID, data.AD$NCFB),
                       group         = rep(data.AD$group, data.AD$NCFB),
                       meanBaseline  = rep(data.AD$MeanBaseline, data.AD$NCFB),
                       sdBaseline    = rep(data.AD$sdBaseline, data.AD$NCFB),
                       meanPost      = rep(data.AD$MeanFU, data.AD$NCFB),
                       sdPost        = rep(data.AD$sdFU, data.AD$NCFB),
                       correlation   = rep(data.AD$Correlation,data.AD$NCFB))

set.seed(123456)
data.IPD$ytmp1 <- rnorm(nrow(data.IPD),0,1)
set.seed(7891011)
data.IPD$ytmp2 <- rnorm(nrow(data.IPD),0,1)

# Standardize ytmp1 and ytmp2, calculate correlation between ytmp1 and ytmp2, 
# and the residuals of regressing ytmp2 on ytmp1
# per study and group

data.IPD2 <- NULL
for(study in unique(data.IPD$study))
{   for (group in unique(data.IPD$group))
{ datatmp     <- data.IPD[data.IPD$study==study & data.IPD$group==group,]
# standardized y1tmp
datatmp$ytmp1 <- (datatmp$ytmp1-mean(datatmp$ytmp1))/sd(datatmp$ytmp1)
# standardized y2tmp
datatmp$ytmp2 <- (datatmp$ytmp2-mean(datatmp$ytmp2))/sd(datatmp$ytmp2)
# correlation between y1tmp and y2tmp
cor.ytmp      <- cor(datatmp$ytmp1, datatmp$ytmp2)
# residuals of regression of ytmp2 on ytmp1
resid         <- residuals(lm(ytmp2 ~ ytmp1 - 1 , data = datatmp))
Resid         <- datatmp$ytmp2 - cor.ytmp*datatmp$ytmp1
# coefficient beta of regression of ytmp2 on ytmp1
#coef         <- coef(lm(ytmp2 ~ ytmp1 - 1 , data = datatmp))
data.IPD2     <- rbind( data.IPD2, data.frame(datatmp,cor.ytmp,resid,Resid))
}  
} 

# temporary variable needed to generate the pseudo baseline and pseudo follow-up outcomes
data.IPD2$ytmp3 <- data.IPD2$ytmp1*data.IPD2$correlation + sqrt(1-data.IPD2$correlation^2)*data.IPD2$resid/sqrt(1-data.IPD2$cor.ytmp^2)
# generate pseudo baseline and pseudo follow-up outcomes
data.IPD2$y1    <- data.IPD2$ytmp1*data.IPD2$sdBaseline + data.IPD2$meanBaseline
data.IPD2$y2    <- data.IPD2$ytmp3*data.IPD2$sdPost + data.IPD2$meanPost

# make new dataset, with only relevant variables
data.pseudoIPD <- data.IPD2[,c("study", "group", "y1", "y2")]
#View(data.pseudoIPD) # final pseudo IPD dataset 
rm(data.IPD2, data.IPD)

# Check the mean and sd of y1 and y2, and correlation y1, y2
check <- cbind(aggregate(y1~group+study, data=data.pseudoIPD, mean), 
              aggregate(y2~group+study, data=data.pseudoIPD, mean)[3],
              aggregate(y1~group+study, data=data.pseudoIPD, sd)[3],
              aggregate(y2~group+study, data=data.pseudoIPD, sd)[3],
              as.vector(cbind(by(data.pseudoIPD, data.pseudoIPD[,c("group","study")], function(x) {cor(x$y1,x$y2)}))))

colnames(check) <- c(colnames(check)[1:2], "meany1", "meany2","sdy1", "sdy2","cory1y2")
check
rm(check)

# Pre-step to calculate centered baseline values by study
data.pseudoIPD$meany1bystudy <- ave(data.pseudoIPD$y1, data.pseudoIPD$study)
data.pseudoIPD$y1center      <- data.pseudoIPD$y1 - data.pseudoIPD$meany1bystudy
data.pseudoIPD$groupcenter   <- data.pseudoIPD$group - 0.5
data.pseudoIPD$arm           <- 1000*data.pseudoIPD$study + data.pseudoIPD$group

#----------------------------------------------------------------------------------------------
#                       One-stage pseusdo IPD models 
#----------------------------------------------------------------------------------------------
ctrl <- lmeControl(opt="optim", msMaxIter=100)

#-----------------------------------------------------------------------------------------------
# Study stratified intercept and random treatment effect ANCOVA main effect

# arm and study specific variances estimated  
FRstudyarm    <- lme(y2 ~ y1center + group + as.factor(study) + y1center*as.factor(study), random= ~ -1 + groupcenter|study, weights =varIdent(form=~study|arm), 
                     control=ctrl, data=data.pseudoIPD, method='REML')

# study specific variances estimated 
FRstudy       <-  lme(y2 ~ y1center+ group + as.factor(study) + y1center*as.factor(study) , random= ~ -1 + groupcenter|study, weights =varIdent(form=~1|study), 
                      control=ctrl, data=data.pseudoIPD, method='REML')

# gruop specific variance estimated 
FRgroup       <-  lme(y2 ~ y1center + group+ as.factor(study) + y1center*as.factor(study) , random= ~ -1 + groupcenter|study, weights =varIdent(form=~1|group),
                      control=ctrl, data=data.pseudoIPD, method='REML')

# one residual variance estimated
FRone         <-  lme(y2 ~ y1center + group + as.factor(study) + y1center*as.factor(study) , random= ~-1 + groupcenter|study, control=ctrl, 
                      data=data.pseudoIPD, method='REML')


# Function to collect the results per model calculating Wald-type CIs
groupeffect <- function(results)
{groupeff   <- summary(results)$tTable["group",]
# Add 95% CI extracted fom nlme using t-distribution
CI          <- intervals(results, which="fixed")
df          <- data.frame(CI$fixed)
names       <- c(names(groupeff),"95% CI lwb", "95% CI upb") 
groupeff    <- c(groupeff, df["group",]$lower, df["group",]$upper)
names(groupeff) <- names  
print(groupeff)
# Estimated tau2  
tau2 <-  VarCorr((results)) [1,]
print("Tau2 is")
print(tau2)
}

groupeffect(FRstudyarm)
groupeffect(FRstudy)
groupeffect(FRgroup)
groupeffect(FRone)

#-----------------------------------------------------------------------------------------------
# Study stratified intercept and random treatment effect ANCOVA interaction effect

# arm and study specific variances estimated
FRstudyarmInt <- lme(y2 ~ y1center*as.factor(study) + y1center*group + group:meany1bystudy, random= ~ -1 + groupcenter|study, weights =varIdent(form=~study|arm), 
                     control=ctrl, data=data.pseudoIPD, method='REML')

# study specific variances estimated  
FRstudyInt    <-   lme(y2 ~ y1center*as.factor(study) + y1center*group + group:meany1bystudy, random= ~ -1 + groupcenter|study, weights =varIdent(form=~1|study), 
                      control=ctrl, data=data.pseudoIPD, method='REML')

# group specific variances estimate
FRgroupInt    <-   lme(y2 ~ y1center*as.factor(study) + y1center*group + group:meany1bystudy, random= ~ -1 + groupcenter|study, weights =varIdent(form=~1|group), 
                       control=ctrl, data=data.pseudoIPD, method='REML')

# one residual variance estimated
FRoneInt      <-   lme(y2 ~ y1center*as.factor(study) + y1center*group + group:meany1bystudy , random= ~ -1 + groupcenter|study, control=ctrl, 
                       data=data.pseudoIPD, method='REML')

# Function to collect the within-trial interactions results per model using Wald-type CIs
within_trial <- function(results)
{inteff <- summary(results)$tTable["y1center:group",]
# Add 95% CI extracted fom nlme using t-distribution
CI        <- intervals(results, which="fixed")
df        <- data.frame(CI$fixed)
names     <- c(names(inteff),"95% CI lwb", " 95% CI upb") 
inteff    <- c(inteff, df["y1center:group",]$lower, df["y1center:group",]$upper)
names(inteff) <- names  
print(inteff)
# Estimated tau2  
tau2 <-  VarCorr((results)) [1,]
print("Tau2 is")
print(tau2)
}

within_trial(FRstudyarmInt)
within_trial(FRstudyInt)
within_trial(FRgroupInt)
within_trial(FRoneInt)

# Function to collect the across-trial interactions results per model using Wald-type CIs
across_trial <- function(results)
{inteff <- summary(results)$tTable["group:meany1bystudy",]
# Add 95% CI extracted fom nlme using t-distribution
CI        <- intervals(results, which="fixed")
df        <- data.frame(CI$fixed)
names     <- c(names(inteff),"95% CI lwb", " 95% CI upb") 
inteff    <- c(inteff, df["group:meany1bystudy",]$lower, df["group:meany1bystudy",]$upper)
names(inteff) <- names  
print(inteff)
# Estimated tau2  
tau2 <-  VarCorr((results)) [1,]
print("Tau2 is")
print(tau2)
}

across_trial(FRstudyarmInt)
across_trial(FRstudyInt)
across_trial(FRgroupInt)
across_trial(FRoneInt)

#----------------------------------------------------------------------------------------------
#                       Two-stage pseusdo IPD models 
#----------------------------------------------------------------------------------------------

# ANCOVA per study on pseudo IPD for subsequent two-stage MA

coef_ancova <- NULL
se_ancova   <- NULL

for (i in unique(data.pseudoIPD$study ))
{         fit <- lm(y2~ y1 + group, data.pseudoIPD[data.pseudoIPD$study==i,])
coef_ancova   <- rbind(coef_ancova,fit$coefficients) 
se_ancova     <- rbind(se_ancova,sqrt(diag(vcov(fit))))
}

# Prepare data for two stage MA
two_stageMA <- data.frame(study=unique(data.pseudoIPD$study), coef_group=coef_ancova[,"group"],
                          secoef_group = se_ancova[,"group"])

# Run aggregate meta-analysis 
MA <- rma(yi <- coef_group, sei=secoef_group, slab=study, method="REML", data=two_stageMA, knha=TRUE)
summary(MA); forest(MA)

#----------------------------------------------------------------------------------------------
# Ancovas per study with interaction of baseline and treatment effect on pseudo IPD for subsequent two-stage MA

coef_ancova_int <- NULL
se_ancova_int   <- NULL

for (i in unique(data.pseudoIPD$study ))
{         fit     <- lm(y2~ y1center + group + y1center*group, data.pseudoIPD[data.pseudoIPD$study==i,])
coef_ancova_int   <- rbind(coef_ancova_int,fit$coefficients) 
se_ancova_int     <- rbind(se_ancova_int, sqrt(diag(vcov(fit))))
}

# Prepare data for two stage MA
two_stageMA_int <- data.frame(study=unique(data.pseudoIPD$study), coef_group=coef_ancova_int[,"y1center:group"],
                              secoef_group = se_ancova_int[,"y1center:group"])

# Run aggregate meta-analysis 
MA_int <- rma(yi <- coef_group,sei=secoef_group, slab=study, method="REML", data=two_stageMA_int)
summary(MA_int); forest(MA_int)
