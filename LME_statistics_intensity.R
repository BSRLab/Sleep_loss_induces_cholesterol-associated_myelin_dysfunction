#R script for data analysis
################################################################

# Load the package for fitting linear mixed effect models - may need to be installed
require(lme4)
library(lmerTest)
library(readxl)
library(tidyverse)
library(sjPlot)
library(jtools)
library(emmeans)
# Load the package for multiple comparisons - may need to be installed
require(multcomp)
require(Rcpp)
library(psych)

#set the directory to upload dataset

df = read_excel('GFAP_BODIPY.xlsx')


#######################################
#make Condition a factor and sleep the reference level for Condition
df$Condition<- as.factor(df$Condition)
df$Rat<- as.factor(df$Rat)

df <- within(df, Condition <- relevel(Condition, ref = '0'))

#check normality of raw data. 
qqnorm(df$Intensity,main="Q-Q Plot Intensity",cex.main=1)
qqline(df$Intensity)
hist(df$Intensity)
shapiro.test(df$Intensity)

# if necessary Use logarithm transformation of the raw data 
#df$Intensity = log(df$Intensity)

# or square root transformation
#df$Intensity = sqrt(df$Intensity)

#check normality of transformed data
#qqnorm(df$Intensity)
#qqline(df$Intensity)
#hist(df$Intensity)
df$Intensity


# 1st step: compare model with only Condition
m1 = lmer(Intensity~1+(1|Rat), data=df, na.action=na.omit, REML=FALSE)
m2 = lmer(Intensity~Condition + (1|Rat), data=df, na.action=na.omit, REML=FALSE)
anova(m1, m2)

#check the normality of residuals of resulting models based on scatterplot and q-q plot of estimated residuals
plot(residuals(m2),main="Residual plot")
qqnorm(residuals(m2),main="Q-Q plot for Residuals")
qqline(residuals(m2))
plot(fitted(m2), residuals(m2))
hist(residuals(m2))
shapiro.test(residuals(m2))
