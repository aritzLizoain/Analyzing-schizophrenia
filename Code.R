# --------------------------------------------------------------------------
# LOAD LIBRARIES AND DATASET
# --------------------------------------------------------------------------

# Libraries
library(data.table)
library(dplyr) 
library(lme4) 
library(nlme)
library(dplyr) 
library("DHARMa")

# Dataset
directory = "" # add file path
setwd(directory)
data = read.table("schizo.txt", sep=" ", header=T)

# --------------------------------------------------------------------------
# EXPLORATORY DATA ANALYSIS
# --------------------------------------------------------------------------

# Filter: keep only individuals with information for all 12 months
data = setDT(data)[, if(.N==12) .SD, by = ID]

# Analyze the general evolution of thought disorders
# Define colors
gray = '#838383'
red = '#ffcccb'
green = '#A6FFAB'
# Plot
boxplot(data$MONTH~data$Y, horizontal=T, col=c(green,red), main='', xlab=expression(bold('Month')), ylab="Thought disorders", names=c(expression(bold('Absence')), expression(bold("Presence"))))
# Add text to the plot
text(x = boxplot.stats(data$MONTH[data$Y=='0'])$stats, labels = boxplot.stats(data$MONTH[data$Y=='0'])$stats, y = 0.5, col=c(gray, gray, 'red', gray, gray))
text(x = boxplot.stats(data$MONTH[data$Y=='1'])$stats, labels = boxplot.stats(data$MONTH[data$Y=='1'])$stats, y = 2.5, col=c(gray, gray, 'red', gray, gray))
text(x = round(mean(na.omit(data$MONTH[data$Y=='1'])),2), labels = paste("Mean:", round(mean(na.omit(data$MONTH[data$Y=='1'])),2)), y = 1.51, col=gray, font = 4)
text(x = round(mean(na.omit(data$MONTH[data$Y=='0'])),2), labels = paste("Mean:", round(mean(na.omit(data$MONTH[data$Y=='0'])),2)), y = 1.49, col=gray, font = 4)
text(x = 8.5, labels = paste("Correlation:", round(cor(data$Y, data$MONTH),2)), y = 1.8 , col='red', font = 4)

# Outlier detection
# Create necessary variables
beginning = numeric(length(unique(data$ID)))
end = beginning
counter = 1
half_length = length(data$Y[data$ID==data$ID[1]])/2
# Loop for all patients
for (id in unique(data$ID)){
  # Number of months with thought disorders in the first 6 months
  beginning[counter] = sum(data$Y[data$ID==id][1:half_length])
  # Number of months with thought disorders in the last 6 months
  end[counter] = sum(data$Y[data$ID==id][(half_length+1):(2*half_length)])
  counter = counter + 1
}
# Color potential outliers in red
colors = rep('black', length(unique(data$ID)))
colors[end>beginning]='red'
# Plot (jitter for visualization purposes)
plot(x=jitter(beginning,.75), y=jitter(end,.75), col=colors, pch=20, cex=2, xlab=expression(bold('Months 0-5')), ylab=expression(bold('Months 6-11')), main='Number of months with thought disorders')
# Add line to the plot
abline(a=0,b=1, col="red", lwd=3) # It is expected a lower number in the last 6 months

# Check the detected potential outliers
unique(data$ID)[end>beginning] # (ID=23 and ID=50))
data[data$ID==23,]
data[data$ID==50,]

# Remove outliers (if considered necessary)
# data = filter(data, ID != 23 & ID!=50)

# Analyze the influence of age and gender
# Function to plot the monthly means
plot_means = function(data, lty){
  table = table(data$MONTH, data$Y)
  means = table[,2]/length(unique(data$ID))
  lines(x=0:11, means, lwd=3, lty=lty)
  # text(x=0:11, labels = round(means,2), col=col, y = means+0.05, font = 2) # Add values to plot
}

# GENDER
# Proportion of patients having no thought disorders as background
cdplot(as.factor(Y) ~ MONTH, data=data, col=c(green,red), main="Evolution depending on GENDER", xlab=expression(bold("Month")), xlim=c(-0.35,11.35), ylab="Thought disorders", yaxlabels = c(expression(bold('Absence')), expression(bold("Presence"))))
# Subset with only male patients
only_male = data[data$GENDER=='male',]
# Subset with only female patients
only_female = data[data$GENDER=='female',]
# Add monthly means for male patients
plot_means(data=only_male, lty=1)
# Add monthly means for female patients
plot_means(data=only_female, lty=2)
# Add legend
legend(x=7, y=0.65, legend=c("Male", "Female"), lty=1:2, cex=1.25, box.lty=0, lwd=2)

# AGE
# Proportion of patients having no thought disorders as background
cdplot(as.factor(Y) ~ MONTH, data=data, col=c(green, red), main="Evolution depending on AGE", xlab=expression(bold("Month")), xlim=c(-0.35,11.35), ylab="Thought disorders", yaxlabels = c(expression(bold('Absence')), expression(bold("Presence"))))
# Subset with only young patients
only_young = data[data$AGE=='young',]
# Subset with only old patients
only_old = data[data$AGE=='old',]
# Add monthly means for young patients
plot_means(data=only_young, lty=1)
# Add monthly means for old patients
plot_means(data=only_old, lty=2)
# Add legend
legend(x=7, y=0.65, legend=c("Young", "Old"), lty=1:2, cex=1.25, box.lty=0, lwd=2)

# --------------------------------------------------------------------------
# MODEL FORMULATION, MODEL SELECTION, AND PARAMETER ESTIMATION
# --------------------------------------------------------------------------

# GLMM 
# First model
glmm1 = glmer(Y ~ MONTH + AGE + GENDER + I(dummy(AGE)*MONTH) + I(dummy(GENDER)*MONTH) + (1|ID), family=binomial(link="logit"), data=data, nAGQ=14)

# Add interaction between AGE and GENDER
glmm2 =  update(glmm1, . ~ . + I(dummy(AGE)*dummy(GENDER)))

# Simplify by removing AGE and GENDER
glmm3 =  update(glmm1, . ~ . - AGE - GENDER)

# Compare
anova(glmm1, glmm2, glmm3)

# Test significance of polynomial components
degree = 3 # linear + quadratic + cubic
glmm3_poly =  update(glmm3, . ~ . + poly(I(as.numeric(GENDER)*MONTH), degree) + poly(I(as.numeric(AGE)*MONTH), degree), nAGQ=0)
# Warning appears because we already had the degree 1 polynomial, so two columns are dropped (one for age and one for gender)
# nAGQ=0 in order to avoid convergence failure

# Compare
anova(glmm3, glmm3_poly)

# Define selected model
final_model = glmm3

# ---------------------------------------------------------------------------------------------------
# MODEL VALIDATION
# ---------------------------------------------------------------------------------------------------

# Assumption: random effects come from a normal distribution, with mean 0, and are uncorrelated.
# Get random effects
random_effects = random.effects(final_model)$ID$`(Intercept)`

# Color in red the largest five random effects
n = length(random_effects)
colors = rep('black', n)
for (i in 1:n){
  if (random_effects[i]==sort(random_effects,partial=n)[n] | random_effects[i]==sort(random_effects,partial=n-1)[n-1] | random_effects[i]==sort(random_effects,partial=n-2)[n-2] | random_effects[i]==sort(random_effects,partial=n-3)[n-3] | random_effects[i]==sort(random_effects,partial=n-4)[n-4]){
    colors[i] = 'red'
  }
}

# Plot
plot(random_effects~unique(data$ID), col=colors, pch=20, cex=2, xlab=expression(bold('Patient ID')), ylab=expression(bold('Random effects')), main='')
# Add mean value line
abline(a=mean(random_effects),b=0, col="red", lwd=3)
par(mfrow=c(1,2))
# ACF
acf(random_effects, main='', xlab=expression(bold('Lag')), ylab=expression(bold('ACF')), lwd=2)
# QQ-plot
qqnorm(random_effects, main="", xlab=expression(bold('Theoretical Quantiles')), ylab=expression(bold('Sample Quantiles')), lwd=2, pch=16)
qqline(random_effects, lwd=2, col='red')
par(mfrow=c(1,1))
# hist(random_effects)

# Check the individuals with the largest five random effects (positive values)
data[data$ID==unique(data$ID)[random_effects == sort(random_effects,partial=n)[n]],] # 68
data[data$ID==unique(data$ID)[random_effects == sort(random_effects,partial=n-1)[n-1]],] # 27
data[data$ID==unique(data$ID)[random_effects == sort(random_effects,partial=n-2)[n-2]],] # 14
data[data$ID==unique(data$ID)[random_effects == sort(random_effects,partial=n-3)[n-3]],] # 33
data[data$ID==unique(data$ID)[random_effects == sort(random_effects,partial=n-4)[n-4]],] # 42

# DHARMa residual analysis
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

# Dispersion problems
# testDispersion(final_model)

# Simulate residuals
simulationOutput = simulateResiduals(fittedModel = final_model, plot = T)

# New model without the patients with the largest five random effects
data_new = filter(data, ID != 68 & ID != 27 & ID != 14 & ID != 33 & ID != 42)
glmm3_new = glmer(Y ~ MONTH + I(dummy(AGE)*MONTH) + I(dummy(GENDER)*MONTH) + (1|ID), family=binomial(link="logit"), data=data_new, nAGQ=14)
# DHARMa residual analysis for the new model
simulationOutput_new = simulateResiduals(fittedModel = glmm3_new, plot = T)

# ---------------------------------------------------------------------------------------------------
# PREDICTIONS
# ---------------------------------------------------------------------------------------------------

# Function to obtain log probability, probability and prediction for a new observation
predict = function(random_intercept, gender, age, month){
  # Get model coefficients
  coeffs = summary(final_model)$coefficients[,1]
  # Define covariates of model
  covariates = c(1,month,age*month,gender*month)
  # Get probability of thought disorders
  log_P = sum(coeffs*covariates)+random_intercept
  P = exp(log_P)/(1+exp(log_P))
  # Set threshold to make a decision
  threshold = 0.5
  prediction = ifelse(P>threshold,1,0)
  # Show prediction
  # print(list(Log_Probability = log_P, Probability=P, Prediction=prediction))
  return(list(Log_Probability = log_P, Probability=P, Prediction=prediction))
}

# Define patient data
random_intercept = rnorm(1, mean=mean(random_effects), sd=sd(random_effects))
gender = 1 # 0='female', 1='male'
age = 1 # 0='old', 1='young'
month = round(runif(1, min=0, max=11),0)

# Predict evolution of one month
male_evolution = predict(random_intercept, 1, age, month+1)$Log_Probability - predict(random_intercept, 1, age, month)$Log_Probability
female_evolution = predict(random_intercept, 0, age, month+1)$Log_Probability - predict(random_intercept, 0, age, month)$Log_Probability
male_evolution - female_evolution # equal to the coefficient for GENDER*MONTH
