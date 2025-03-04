library(car)
library(pROC)
library(MASS)
library(class)
library(boot)
library(leaps)
library(glmnet)
library(ISLR2)
library(corrplot)

# read.csv(file = "~/koi.csv")

# Check for duplicates
duplicates <- sum(duplicated(koi))
duplicates

# We have no duplicates, that means that every row corresponds to a different KOI.

# Check for duplicates in the 'kepid' column
duplicates_kepid <- sum(duplicated(koi$kepid))
duplicates_kepid

# But we have 1350 duplicates in kepid. That means that 1350 stars have at least
# 2 KOI associated to them.

# Create a table of frequencies for kepid values
kepid_table <- table(koi$kepid)

# Plot the frequency of kepid values
barplot(kepid_table,
        main = "Frequency of kepid values",
        xlab = "kepid",
        ylab = "Frequency",
        las = 2,  
        cex.names = 0.6)

# By plotting the frequencies of kepid values, we find out that stars have at
# least one KOI associated to them, up to 7. Now we want to find out how common
# for a star is to have n KOI associated to them.

# Count how many times each frequency occurs
frequency_count <- table(kepid_table)

# Print the frequency counts for 1 to 7 times
for (i in 1:7) {
  count <- ifelse(i %in% names(frequency_count), frequency_count[as.character(i)], 0)
  print(paste("Number of kepid values present", i, "times:", count))
}

# Count the number of unique values in the kepid column
unique_kepid_count <- length(unique(koi$kepid))
unique_kepid_count

# This tells us that we have a total of 8214 stars in our dataset.


# Create a table of frequencies for koi_disposition values
disposition_frequency <- table(koi$koi_disposition)

# Plot the distribution of koi_disposition values using a bar plot
barplot(disposition_frequency,
        main = "Distribution of koi_disposition values",
        xlab = "koi_disposition",
        ylab = "Frequency",
        col = c("lightblue", "aquamarine", "bisque"))

# Check for NA values in the dataset
na_count <- sapply(koi, function(x) sum(is.na(x)))
na_count

# We realized that there are no NA values in kepler_name, but by looking at the
# dataset we can see that there are missing values.

# We substitute empty strings with NA values in kepler_name
koi$kepler_name[koi$kepler_name == ""] <- NA
na_count <- sapply(koi, function(x) sum(is.na(x)))
na_count

# Now we can see that we have actually 6819 NAs in kepler_name.
# What about other columns?

# Define a vector of ghost values to be replaced with NA
ghost_values <- c("", "NULL", "-", " ", "N/A")  

# Function to replace ghost values with NA in a column
replace_ghost_values <- function(column) {
  column[column %in% ghost_values] <- NA
  return(column)
}

# Apply the function to each column
koi <- as.data.frame(lapply(koi, replace_ghost_values))

# Verify that ghost values have been replaced with NA
na_count <- sapply(koi, function(x) sum(is.na(x)))
na_count

# We found out that there were no other ghost NAs in the dataset.


# The majority of NAs are in kepler_name. That is to be expected, as KOI only
# get named after they have been confirmed.

# Count the number of CONFIRMED values in the koi_disposition column
confirmed_count <- sum(koi$koi_disposition == "CONFIRMED")
confirmed_count

# We have 2743 CONFIRMED KOIs, 9564 (total n of KOIs) - 6819 (KOIs without a name)
# = 2745

# We find in this way that there are 2 KOIs, either CANDIDATE or
# FALSE POSITIVE, that have a name.

# Filter the rows where koi_disposition is either "CANDIDATE" or "FALSE POSITIVE" and kepler_name is not NA
filtered_koi <- subset(koi, (koi_disposition %in% c("CANDIDATE", "FALSE POSITIVE")) & !is.na(kepler_name))

# Retrieve the kepler_name column for these rows
kepler_names <- filtered_koi$kepler_name
kepler_names

# Looking at the dataset, we find that "Kepler-469 b" and "Kepler-503 b"
# are both FALSE POSITIVE.

# Let's continue with our data-cleaning. There are 1510 NAs in koi_score.

# The koi_score variable represents a quantitative measure of the likelihood
# that a Kepler Object of Interest (KOI) is a true exoplanet.
# It takes values ranging between 0 and 1.

# Now, we can see that some variables have more or less the same amount of NAs.
# A good idea would be to eliminate rows that present NAs in one variable, 
# and then check if also other NAs disappear.

# Remove rows with NA in the koi_impact column
koi <- koi[!is.na(koi$koi_impact), ]

na_count <- sapply(koi, function(x) sum(is.na(x)))
na_count


# We are left with 255 NAs in both koi_tce_plnt_num and koi_tce_delivname,
# and one NA in koi_kepmag.

# Both koi_tce_plnt_num and koi_tce_delivname are not useful variables
# for our analysis, since they are used to classify the Threshold Crossing Event.

# koi_kepmag is a useful variable, which represents the Kepler magnitude
# of the host star. Let's remove the associated NA value.

# Remove rows with NA in the koi_impact column
koi <- koi[!is.na(koi$koi_kepmag), ]

na_count <- sapply(koi, function(x) sum(is.na(x)))
na_count


# Separate rows with NA in koi_score into a new dataset
koi_na <- koi[is.na(koi$koi_score), ]

# Remove rows with NA in koi_score from the original dataset
koi <- koi[!is.na(koi$koi_score), ]

na_count <- sapply(koi, function(x) sum(is.na(x)))
na_count


# We are left with 5267 NAs in kepler_name, and a total of 7994 objects.
# There are no other NAs in the dataset.




### EDA


## koi_disposition
# Create a table of frequencies for koi_disposition values
disposition_frequency <- table(koi$koi_disposition)

# Plot the distribution of koi_disposition values using a bar plot
barplot(disposition_frequency,
        main = "Distribution of koi_disposition values",
        xlab = "koi_disposition",
        ylab = "Frequency",
        col = c("lightblue", "aquamarine", "bisque"))

## koi_pdisposition
# Create a table of frequencies for koi_pdisposition values
pdisposition_frequency <- table(koi$koi_pdisposition)

# Plot the distribution of koi_disposition values using a bar plot
barplot(pdisposition_frequency,
        main = "Distribution of koi_disposition values",
        xlab = "koi_disposition",
        ylab = "Frequency",
        col = c("lightblue", "bisque"))

## koi_score
# Plot the distribution of koi_score using a histogram
hist(koi$koi_score,
     main = "Distribution of koi_score",
     xlab = "koi_score",
     ylab = "Frequency",
     col = "lightblue",  
     breaks = 30)


## koi_fpflag_nt
# Create a table of frequencies for koi_fpflag_nt values
fpflag_nt_frequency <- table(koi$koi_fpflag_nt)

# Plot the distribution of koi_fpflag_nt values using a bar plot
barplot(fpflag_nt_frequency,
        main = "Distribution of koi_fpflag_nt values",
        xlab = "koi_fpflag_nt",
        ylab = "Frequency",
        col = c("lightblue", "aquamarine", "bisque"),  
        names.arg = c("0", "1", "465"))

# The koi_fpflag_nt variable is a binary flag that indicates whether the
# Kepler Object of Interest (KOI) was identified as a false positive due to
# not having a transit-like signal. This helps in distinguishing between
# true planetary transits and other types of signals that do not resemble
# planetary transits.

# We have an error in our dataset. Since the koi_disposition is CONFIRMED,
# we know that the right value would have been 0.

# Change the value 465 to 0 in the koi_fpflag_nt column
koi$koi_fpflag_nt[koi$koi_fpflag_nt == 465] <- 0


## koi_fpflag_ss
# Create a table of frequencies for koi_fpflag_nt values
fpflag_ss_frequency <- table(koi$koi_fpflag_ss)

# Plot the distribution of koi_fpflag_nt values using a bar plot
barplot(fpflag_ss_frequency,
        main = "Distribution of koi_fpflag_ss values",
        xlab = "koi_fpflag_ss",
        ylab = "Frequency",
        col = c("lightblue", "aquamarine"),  
        names.arg = c("0", "1"))


## koi_fpflag_co
# Create a table of frequencies for koi_fpflag_co values
fpflag_co_frequency <- table(koi$koi_fpflag_co)

# Plot the distribution of koi_fpflag_nt values using a bar plot
barplot(fpflag_co_frequency,
        main = "Distribution of koi_fpflag_co values",
        xlab = "koi_fpflag_co",
        ylab = "Frequency",
        col = c("lightblue", "aquamarine"),  
        names.arg = c("0", "1"))


## koi_fpflag_ec
# Create a table of frequencies for koi_fpflag_co values
fpflag_ec_frequency <- table(koi$koi_fpflag_ec)

# Plot the distribution of koi_fpflag_nt values using a bar plot
barplot(fpflag_ec_frequency,
        main = "Distribution of koi_fpflag_ec values",
        xlab = "koi_fpflag_ec",
        ylab = "Frequency",
        col = c("lightblue", "aquamarine"),  
        names.arg = c("0", "1"))

## koi_period
# Plot the distribution of koi_period using a histogram
hist(koi$koi_period,
     main = "Distribution of koi_period",
     xlab = "koi_period",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

# We can see how koi_period has a lot of outliers, let's use the boxplot plot
# to better visualize outliers are distributed

# Plot the distribution of koi_period using a boxplot
boxplot(koi$koi_period,
        main = "Boxplot of koi_period",
        ylab = "koi_period",
        col = "lightblue",  # Optional: add color to the box
        horizontal = T)

## koi_time0bk
# The koi_time0bk variable represents the time of the first detected transit
# center for a Kepler Object of Interest (KOI). It's a Barycentric Kepler
# Julian Date, it is used for classifying the periodic dimming of KOIs,
# but it is not useful for our analysis.

## koi_impact
hist(koi$koi_impact,
     main = "Distribution of koi_impact",
     xlab = "koi_impact",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

boxplot(koi$koi_impact,
        main = "Boxplot of koi_impact",
        ylab = "koi_impact",
        col = "lightblue",  
        horizontal = TRUE)

## koi_duration
hist(koi$koi_duration,
     main = "Distribution of koi_duration",
     xlab = "koi_duration",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

boxplot(koi$koi_duration,
        main = "Boxplot of koi_duration",
        ylab = "koi_duration",
        col = "lightblue",  
        horizontal = TRUE)

## koi_depth
hist(koi$koi_depth,
     main = "Distribution of koi_depth",
     xlab = "koi_depth",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

boxplot(koi$koi_depth,
        main = "Boxplot of koi_depth",
        ylab = "koi_depth",
        col = "lightblue",  
        horizontal = TRUE)

## koi_prad
hist(koi$koi_prad,
     main = "Distribution of koi_prad",
     xlab = "koi_prad",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 100)

## koi_teq
hist(koi$koi_teq,
     main = "Distribution of koi_teq",
     xlab = "koi_teq",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

## koi_insol
hist(koi$koi_insol,
     main = "Distribution of koi_insol",
     xlab = "koi_insol",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 100)

# koi_model_snr
hist(koi$koi_model_snr,
     main = "Distribution of koi_model_snr",
     xlab = "koi_model_snr",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

# koi_tce_plnt_num and koi_tce_delivname are used to classify the dimness events
# they are not useful for our analysis


## koi_steff
hist(koi$koi_steff,
     main = "Distribution of koi_steff",
     xlab = "koi_steff",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

## koi_slogg
hist(koi$koi_slogg,
     main = "Distribution of koi_slogg",
     xlab = "koi_slogg",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

# koi_steff and koi_slogg presents less outliers, so their distributions 
# is more or less similar to a Normal one

## koi_srad
hist(koi$koi_srad,
     main = "Distribution of koi_srad",
     xlab = "koi_srad",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

# ra and dec are coordinates, they are not useful for our analysis

## koi_kepmag
hist(koi$koi_kepmag,
     main = "Distribution of koi_kepmag",
     xlab = "koi_kepmag",
     ylab = "Frequency",
     col = "lightblue",
     breaks = 30)

# Let's remove variables that are not of interest.
koi$ra <- NULL
koi$dec <- NULL
koi$koi_tce_plnt_num <- NULL
koi$koi_tce_delivname <- NULL
koi$koi_time0bk <- NULL

# Lots of variables are connected to each other through formulas,
# and also by their nature. 


cov_matrix <- round(var(koi[,-1:-5]), 2)
cov_matrix

corr_matrix <- round(cor(koi[,-1:-5]), 3)

corrplot(corr_matrix, method = 'circle', type = "upper",
         tl.col = "black", cl.pos = "n")

# Some variables show significant (negative or positive) linear relationship
# Let's explore them.

plot(koi$koi_insol, koi$koi_srad, pch=20)
plot(koi$koi_slogg, koi$koi_kepmag, pch=20)
plot(koi$koi_teq, koi$koi_slogg, pch=20)
plot(koi$koi_slogg, koi$koi_srad, pch=20)


plot(koi$koi_fpflag_co, koi$koi_fpflag_ec, pch=20)
correlation <- cor(koi$koi_fpflag_co, koi$koi_fpflag_ec)
correlation




koi.r <- koi


koi.r$kepid <- NULL
koi.r$kepoi_name <- NULL
koi.r$kepler_name <- NULL
koi.r$koi_disposition <- NULL
koi.r$koi_pdisposition <- NULL

write.csv(koi.r, "koi_r.csv")


#------------

# Multiple Linear Regression

mlr.out <- lm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
              + koi_fpflag_ec + koi_period + koi_impact + koi_duration
              + koi_depth + koi_prad + koi_teq + koi_insol
              + koi_model_snr + koi_steff + koi_slogg + koi_srad
              + koi_kepmag, data=koi.r)
summary(mlr.out)

vif_mlr.out <- vif(mlr.out)
vif_mlr.out

# Residuals vs. Fitted Plot
plot(mlr.out, which = 1)

# Normal Q-Q Plot
plot(mlr.out, which = 2)

# Scale-Location Plot
plot(mlr.out, which = 3)

# Residuals vs. Leverage Plot
plot(mlr.out, which = 5)

# The variables koi_srad and koi_kepmag are not significant.
# The model explains about 71.29% of the variance in koi_score.
# Adjusted R-squared is slightly lower than the multiple R-squared,
# indicating that the model is well-fitted.
# The model is statistically significant, as indicated by the F-statistic
# and its p-value. The coefficient estimate are very small.

# Removing variables with VIF > 2.4
mlr.out.R <- lm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
              + koi_fpflag_ec + koi_period + koi_impact + koi_duration
              + koi_depth + koi_prad + koi_insol
              + koi_model_snr + koi_steff
              + koi_kepmag, data=koi.r)
summary(mlr.out.R)

# ANOVA
anova(mlr.out.R, mlr.out)
# The RSS tells us that the excluded variable do not improve the fit of the model
# that much, at the same time, the extremely small p-value suggests we should not
# exclude the 3 variables from the model.


# Residuals vs. Fitted Plot
plot(mlr.out.R, which = 1)

# Normal Q-Q Plot
plot(mlr.out.R, which = 2)

# Scale-Location Plot
plot(mlr.out.R, which = 3)

# Residuals vs. Leverage Plot
plot(mlr.out.R, which = 5)

# Removing the 3 variables did not improve the plots. The problems are probably
# due to heteroschedasticity, and from the fact that the relationship between
# the variables and koi_score is not linear (and maybe, is not present at all).



mlr.out.nb <- lm(koi_score ~ koi_period + koi_impact + koi_duration
              + koi_depth + koi_prad + koi_teq + koi_insol
              + koi_model_snr + koi_steff + koi_slogg + koi_srad
              + koi_kepmag, data=koi.r)
summary(mlr.out.nb)

vif_mlr.out.nb <- vif(mlr.out.nb)
vif_mlr.out.nb

# Residuals vs. Fitted Plot
plot(mlr.out.nb, which = 1)

# Normal Q-Q Plot
plot(mlr.out.nb, which = 2)

# Scale-Location Plot
plot(mlr.out.nb, which = 3)

# Residuals vs. Leverage Plot
plot(mlr.out.nb, which = 5)


# R-squared dropped to 0.3208, indicating a moderate fit.
# koi_steff and koi_srad are not significant.
# Overall, the model is statistically significant, as indicated by the very
# low p-value


mlr.out.b <- lm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co 
                 + koi_fpflag_ec, data = koi.r)
summary(mlr.out.b)

# All four predictors are highly significant, indicating a strong negative
# relationship with koi_score. The model explains about 67.13% of the variance
# in koi_score, indicating a good fit

# Residuals vs. Fitted Plot
plot(mlr.out.b, which = 1)
# The residuals do not seem to be randomly scattered around zero. Instead,
#there is a clear pattern, with residuals first decreasing and then increasing.
# This indicates that the relationship between the predictors and the response
# variable might not be linear.
# Also, the varying spread of the residuals suggests heteroscedasticity.
# This means that the variance of the errors is not constant, which violates
# one of the assumptions of linear regression.

# Normal Q-Q Plot
plot(mlr.out.b, which = 2)
#The points form an S-shaped curve, which indicates that the residuals are not
# normally distributed. There is a clear deviation from normality in both tails,
# and the middle part shows a pattern as well.

# Scale-Location Plot
plot(mlr.out.b, which = 3)
# The red line shows an increasing trend with the fitted values. This confirms
# that the variance of the residuals is not constant, suggesting
# heteroscedasticity.

# Residuals vs. Leverage Plot
plot(mlr.out.b, which = 5)
# There are some points with higher leverage values, particularly around 0.0015
# to 0.0025 on the leverage axis. These points have both high leverage and larg
# residuals.


vif_mlr.out.b <- vif(mlr.out.b)
vif_mlr.out.b

# All the VIF values are between 1 and 2, indicating low to moderate
# multicollinearity. This was particularly important to test since fpflag
# variables are binary, and we had to make sure to avoid the dummy variables
# trap. What we can evaluate by looking at the dataset is, when at least one
# variable is flagged, the koi_score drops to 0 or close, but more than one
# variable could be flagged. Also, if no variables are flagged, koi_score
# could still be very low or close to zero. In a sense, fpflag variables are not
# providing overlapping or complete information. They provide information about
# why a planet may have a low score, but they do not provide information
# about how to classify a planet that has not been flagged.


# ANOVA
anova(mlr.out.b, mlr.out)
anova(mlr.out.nb, mlr.out)



# Investigating the relationship between koi_score and the dipendent variables
boxplot(koi_score ~ koi_fpflag_co, data = koi.r, col = "lightblue")

boxplot(koi_score ~ koi_fpflag_nt, data = koi.r, col ="aquamarine")

boxplot(koi_score ~ koi_fpflag_ss, data = koi.r, col = "bisque")

boxplot(koi_score ~ koi_fpflag_ec, data = koi.r, col = "lightcoral")


plot(koi.r$koi_score, koi.r$koi_period)
# koi_period values vary significantly, with many values concentrated near 0
# and a few observations having very high periods (up to over 800 days).
# There doesn't appear to be a clear linear relationship between koi_score
# and koi_period.
# The distribution of koi_period values remains relatively constant across
# different levels of koi_score, indicating that koi_period might not be
# strongly dependent on koi_score.
koi.r$log_koi_period <- log(koi.r$koi_period + 1)
plot(koi.r$koi_score, koi.r$log_koi_period)

koi.r$log_koi_score <- log(koi.r$koi_score)
plot(koi.r$log_koi_score, koi.r$log_koi_period)
plot(koi.r$log_koi_score, koi.r$koi_period)
plot(log(koi.r$koi_score + 1e-8 / 1 - koi.r$koi_score + 1e-8), koi.r$koi_period)

plot(koi.r$koi_score, koi.r$koi_impact)
# koi_impact usually takes values ranging between 0 and 1. When koi_impact is
# higer, it means the planet’s path misses the star's disk, indicating no
# transit would occur.
# koi_impact outliers seem to be associated with lower koi_scores, otherwise
# no clear relationship arises from the graph.
koi.r$log_koi_impact <- log(koi.r$koi_impact + 1)
plot(koi.r$koi_score, koi.r$log_koi_impact)

plot(koi.r$koi_score, koi.r$koi_duration)
# also in this case, no clear relationship arises from the graph,
# but we can see that higher koi_duration leads either to scores of 0 and 1,
# that makes sense since a higer duration (in days) helps in analyzing the
# shape and characteristics of the transit light curve.
# In a sense, higher durations lead to more certain classifications.
koi.r$log_koi_duration <- log(koi.r$koi_duration + 1)
plot(koi.r$koi_score, koi.r$log_koi_duration)

plot(koi.r$koi_score, koi.r$koi_depth)
# koi_depth variable represents the depth of the transit, which is the amount
# by which the star's brightness decreases when the exoplanet passes in front
# of it. This depth is a measure of how much light is blocked by the planet
# and provides insights into the size of the planet relative to the star.

# Transforming depths from ppm into percentages
koi.r$koi_depth_p <- koi.r$koi_depth / 1000
plot(koi.r$koi_score, koi.r$koi_depth_p)

# Applying the log-transform to stabilize the variance
koi.r$log_koi_depth_p <- log(koi.r$koi_depth_p + 1)
plot(koi.r$koi_score, koi.r$log_koi_depth_p)

# Even after both transformation, there still seems to be heteroschedasticity
# but before, it looked like bigger depths were associated strongly with lower
# koi_scores, now it seems there is no relationship between the two variables
# Bigger depths, allow for more certain classification of both types.

plot(koi.r$koi_score, koi.r$koi_prad)
# The koi_prad is a numeric value that indicates the size of the exoplanet
# relative to Earth's radius.
# The graph is difficult to interpret since there is one important outlier that
# flattens other observations.
koi.r$log_koi_prad <- log(koi.r$koi_prad + 1)
plot(koi.r$koi_score, koi.r$log_koi_prad)
# The log transformation of koi_prad has compressed the range of values and
# reduced the impact of extreme outliers. Values appears more uniform across
# different levels of koi_score, indicating that the transformation has
# helped mitigate heteroscedasticity.
# There is still no clear relationship between koi_score and
# log_koi_prad.


plot(koi.r$koi_score, koi.r$koi_teq)
# Convert Kelvin into Celsius
koi.r$celsius_koi_teq <- (koi.r$koi_teq - 273.15)
plot(koi.r$koi_score, koi.r$celsius_koi_teq)
# Also in this case, there seems to be no relationship between koi_score and the
# the estimated equilibrium temperature of the potential planets.
# Higher temperatures lead to more certain classifications, as more luminous star
# are associated with an imorovement in the strength and clarity of the
# transit signal.
koi.r$log_celsius_koi_teq <- log(koi.r$celsius_koi_teq + 1)
plot(koi.r$koi_score, koi.r$log_celsius_koi_teq)

plot(koi.r$koi_score, koi.r$koi_insol)
# koi_insol is a numeric value that indicates the amount of stellar energy the
# exoplanet receives compared to what Earth receives from the Sun.
# Extreme outlier make the graph hard to interpret.
koi.r$log_koi_insol <- log(koi.r$koi_insol + 1)
plot(koi.r$koi_score, koi.r$log_koi_insol)
# Same as before, the amount of stellar energy a planet receives is strictly
# linked with the star's luminosity, the distance from the star and orbital
# eccentricity, variables that also lead to clearer classifications.


plot(koi.r$koi_score, koi.r$koi_model_snr)
# koi_model_snr variable represents the signal-to-noise ratio (SNR) of the
# detected transit signal for a Kepler Object of Interest (KOI).
# This metric quantifies the strength of the transit signal relative to the
# background noise in the observational data.
# We realized that other variables are strictly linked to the SNR, meaning
# we may have a multi-collinearity problem that needs to be addressed.
# Looking at the formula used to calculate the SNR, we see how it is completely
# explained by the transit depth, scaled byt the standard deviation of noise
# in the observational data.
koi.r$log_koi_model_snr <- log(koi.r$koi_model_snr + 1)
plot(koi.r$koi_score, koi.r$log_koi_model_snr)


plot(koi.r$koi_score, koi.r$koi_steff)
# koi_steff variable represents the effective temperature of the host star of a
# Kepler Object of Interest.
# Same considerations as before.
koi.r$celsius_koi_steff <- (koi.r$koi_steff - 273.15)
plot(koi.r$koi_score, koi.r$celsius_koi_steff)

koi.r$log_celsius_koi_steff <- log(koi.r$celsius_koi_steff + 1)
plot(koi.r$koi_score, koi.r$log_celsius_koi_steff)



plot(koi.r$koi_score, koi.r$koi_slogg)
# koi_slogg variable represents the logarithm of the surface gravity of the
# host star of a Kepler Object of Interest. 

plot(koi.r$koi_score, koi.r$koi_srad)
# koi_srad variable represents the radius of the host star of a Kepler Object
# of Interest (KOI), expressed in units of solar radii. The size of the star
#directly influences the interpretation of the transit depth, which is used
# to calculate the exoplanet's radius. Also, combined with other stellar
# parameters, the radius helps estimate the star’s luminosity and effective
# temperature. Same considerations as before.
koi.r$log_koi_srad <- log(koi.r$koi_srad + 1)
plot(koi.r$koi_score, koi.r$log_koi_srad)

plot(koi.r$koi_score, koi.r$koi_kepmag)
# koi_kepmag variable represents the Kepler magnitude of the host star of a
# Kepler Object of Interest (KOI). The Kepler magnitude is a measure of the
# star's brightness as observed by the Kepler space telescope.
# The brightness of the star affects the quality and accuracy of the photometric
# measurements used to detect and characterize exoplanets.
# The idea is that brighter stars are easier to study, but we can see from the
# graph that also dimmer stars (higher kepmag values) lead to more certain
# classifications.


# Let's see if we can improve the fit by transforming the variables:
# linear-log
mlr.linlog.out <- lm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
              + koi_fpflag_ec + log_koi_period + log_koi_impact + log_koi_duration
              + log_koi_depth_p + log_koi_prad + log_celsius_koi_teq + log_koi_insol
              + log_koi_model_snr + log_celsius_koi_steff + koi_slogg + log_koi_srad
              + koi_kepmag, data=koi.r)
summary(mlr.linlog.out)

# log-log
mlr.loglog.out <- lm(log(koi_score+0.1) ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
                  + koi_fpflag_ec + log_koi_period + log_koi_impact + log_koi_duration
                  + log_koi_depth_p + log_koi_prad + log_celsius_koi_teq + log_koi_insol
                  + log_koi_model_snr + log_celsius_koi_steff + koi_slogg + log_koi_srad
                  + koi_kepmag, data=koi.r)
summary(mlr.loglog.out)


# residual standard error is lower for the lin-log model, meaning 
# predictions are closer to the observed values on average.
# At the same time the R-squared is higher for the log-log model, meaning
# the model can explain more variance in the response variable.
# The log-log and complete log-log models (applying log on all variables) are
# equivalent.

# To decide which model is best, let's plot diagnostic plots.
# Residuals vs. Fitted Plot
plot(mlr.linlog.out, which = 1)

# Normal Q-Q Plot
plot(mlr.linlog.out, which = 2)

# Scale-Location Plot
plot(mlr.linlog.out, which = 3)

# Residuals vs. Leverage Plot
plot(mlr.linlog.out, which = 5)


# Residuals vs. Fitted Plot
plot(mlr.loglog.out, which = 1)

# Normal Q-Q Plot
plot(mlr.loglog.out, which = 2)

# Scale-Location Plot
plot(mlr.loglog.out, which = 3)

# Residuals vs. Leverage Plot
plot(mlr.loglog.out, which = 5)

# The graphs are completely equivalent. In order to make it easier to interpret
# the effect on koi_score, we can use the lin-log model.
# Also, there is no improvement from the graph plotted before the log-transform.
# We only solved the outliers/high-leverage points problem.
# That explains the slightly better fit of the models.

# We have tried to fit a multiple linear regression on a complex dataset.
# The relationship between the variables and our koi_score are complicated,
# and a linear model can not explain it fully. Even though our R-squared was
# good enough (around 0.75), most of the dipendent variables is explained
# by the 4 fpflag binary variables.
# This means that the rest of the variables were not fitted
# This could be due to the non-linearity, or heteroschedasticity, coti

# That is why, we decided to fit a logistic regression, not to check if
# but to check if a planet is easy to classifi or not





# Logistic Regression

# To check how certain it is to classify a signal as either a planet or a false positive
koi.r$is_classifiable <- ifelse(koi.r$koi_score < 0.25 | koi.r$koi_score > 0.75,
                                1, 0)
# To classify a signal as a planet or a false positive based on koi_score
koi.r$is_planet <- ifelse(koi.r$koi_score > 0.5, 1, 0)

# To classify a signal as a planet or not based on koi_disposition
koi.r$confirmed <- ifelse(koi.r$koi_disposition == "CONFIRMED", 1, 0)

# To proceed with a logistic regression, we thought about 3 possible response variables
# Our data made us think that our regressors are not strictly linked with
# koi_score, making it hard to predict/classify whether a signal represents a 
# planet, or is a false positive. For middle values of koi_score, our scatterplots
# show how the data is dispersed (seemingly) randomly, but for extreme values of
# koi_score, there seems to be some relation.
# That is why we thought about classifying for easy/difficult classification,
# as a tool to indicate which planets might be smarter to focus on for classification,
# assuming other data will be available for the actual classification.
# But our logistic regression on is_classifiable gave us a neglectable reduction
# from the null to the residual deviance.

# Surprisingly this did not happen with is_planet, which gave us a consistent 
# improvement in reducing the deviance.

# But this variable was imprecise, as it is common in our dataset that for koi_score
# > 0.5 we have lots of candidate or false positive.

# That is why we used koi_disposition to create the response variable.
# Also, the glm gave us a warning message (glm.fit: fitted probabilities numerically 0 or 1 occurred)
# meaning we have complete or quasi-complete separation, that is,
# when some of the predictor variables perfectly separate the response variable into 0s and 1s.
# Also the AIC is larger than it is on is_planet.
# We choose to continue with this model, as it poses for improvements using
# regularization techniques.

confirmed.out <- glm(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
                     + koi_fpflag_ec + koi_period + koi_impact + koi_duration
                     + koi_depth + koi_prad + koi_teq + koi_insol
                     + koi_model_snr + koi_steff + koi_slogg + koi_srad
                     + koi_kepmag, family = binomial, data=koi.r)
summary(confirmed.out)


# Obtain the estimated probabilities 
logistic.prob <- predict(confirmed.out, type="response")

# Assign probabilities > 0.5 to the class 1
logistic.pred <- rep(0, 7994)
logistic.pred[logistic.prob > 0.6] <- 1

conf.matrix <- table(logistic.pred, koi.r$confirmed)
dimnames(conf.matrix) <- list("Predictions" = c(0, 1), "Confirmed" = c(0, 1))
conf.matrix



# Function to compute performance measures
#
# Arguments:
#
# pred.values = vector of predicted values
# true.values = vector of true values
# lab.pos     = label of the positive class
#
perf.measure <- function(true.values, pred.values,  lab.pos = 1){
  #
  # compute the confusion matrix and number of units
  conf.matrix <- table(pred.values, true.values)
  n <- sum(conf.matrix)
  #
  # force the label of positives to be a character string
  lab.pos <- as.character(lab.pos)
  #
  # obtain the label of negatives
  lab <- rownames(conf.matrix)
  lab.neg <- lab[lab != lab.pos]
  #
  # extract relevant quantities from the confusion matrix
  TP <- conf.matrix[lab.pos, lab.pos]
  TN <- conf.matrix[lab.neg, lab.neg]
  FP <- conf.matrix[lab.pos, lab.neg]
  FN <- conf.matrix[lab.neg, lab.pos]
  P     <- TP + FN
  N     <- FP + TN
  P.ast <- TP + FP
  #
  # compute the performance measures
  OER <- (FP+FN)/n
  PPV <- TP/P.ast
  TPR <- TP/P
  F1  <- 2*PPV*TPR/(PPV+TPR)
  TNR <- TN/N
  FPR <- FP/N
  return(list(Overall.Error.Rate = OER, Precision=PPV, Recall=TPR, F1=F1, Specificity=TNR, False.Positive.Rate=FPR))
}

PM <- perf.measure(koi.r$confirmed, logistic.pred,  lab.pos = 1)
PM

# We compared 0.5, 0.66, and 0.6 threshold.
# While the 0.5 threshold gives higer Recall and F1 scores, the 0.66 threshold
# maximizes Precision and Specificity while minimizing False Positives.
# 0.6 threshold offers a baalanced middle-ground between the two previous threshold.



# ROC curve 
roc.out <- roc(koi.r$confirmed, logistic.prob, levels=c(0, 1))

# different ways of plotting the ROC curve
plot(roc.out, print.auc=TRUE) # check values on the x axis


# Threshold that maximizes the sum of Specificity and Sensitivity
coords(roc.out, "best")

# threshold = "best"
best.th <- coords(roc.out, "best")$threshold
logistic.pred.best <- rep(0, 7994)
logistic.pred.best[logistic.prob>best.th] <- 1


confusion.matrix.best <- table(logistic.pred.best, koi.r$confirmed)
dimnames(confusion.matrix.best) <- list("Predictions" = c(0, 1), "Confirmed" = c(0, 1))
confusion.matrix.best

# The best threshold (0.4417261) is more prone to Type I errors, while greatly
# reducing Type II errors.

PM.best <- perf.measure(koi.r$confirmed, logistic.pred.best,  lab.pos = 1)
PM.best




# Selection of predictors

confirmed.out.F <- glm(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
                     + koi_fpflag_ec + koi_period + koi_impact + koi_duration
                     + koi_depth + koi_prad + koi_teq + koi_insol
                     + koi_model_snr + koi_steff + koi_slogg + koi_srad
                     + koi_kepmag, family = binomial, data=koi.r)
summary(confirmed.out.F)

# step 1
confirmed.out.1 <- update(confirmed.out.F, .~.-koi_insol)
summary(confirmed.out.1)

anova(confirmed.out.1, confirmed.out.F, test="Chisq")

# step 2

confirmed.out.2 <- update(confirmed.out.1, .~.-koi_slogg)
summary(confirmed.out.2)

anova(confirmed.out.2, confirmed.out.F, test="Chisq")

# step 3

confirmed.out.3 <- update(confirmed.out.2, .~.-koi_impact)
summary(confirmed.out.3)

anova(confirmed.out.3, confirmed.out.F, test="Chisq")

# step 4

confirmed.out.4 <- update(confirmed.out.3, .~.-koi_duration)
summary(confirmed.out.4)

anova(confirmed.out.4, confirmed.out.F, test="Chisq")

# step 5

confirmed.out.5 <- update(confirmed.out.4, .~.-koi_srad)
summary(confirmed.out.5)

anova(confirmed.out.5, confirmed.out.F, test="Chisq")

# step 6

confirmed.out.6 <- update(confirmed.out.5, .~.-koi_fpflag_ec)
summary(confirmed.out.6)

anova(confirmed.out.6, confirmed.out.F, test="Chisq")

# step 7

confirmed.out.7 <- update(confirmed.out.6, .~.-koi_fpflag_co)
summary(confirmed.out.7)

anova(confirmed.out.7, confirmed.out.F, test="Chisq")

# We removed variables based on their significance.

# Until model 5, we were able to remove variables without increasing much the
# residual deviance, we went from 4898.6 to 4907.8.
# With model 6 we got to 5076.1. With model 7 to 6112.6.

# Even though we were removing variables, the AIC got bigger, especially with
# models 6 and 7.

# The anova F-statistic was always indicating that the reduction in
# deviance was highly statistically significant, starting from model 1.

# We were also not able to resolve the warning message "fitted probabilities
# numerically 0 or 1 occurred".

# Choosing only on the differences in residual deviances, the best model is
# probably the sixth one.


# Training and Validation Set

set.seed(42)

train_size <- floor(0.7 * nrow(koi.r))

train_indices <- sample(seq_len(nrow(koi.r)), size = train_size)

# Split the data into training and validation sets
train_set <- koi.r[train_indices, ]
validation_set <- koi.r[-train_indices, ]


# Classification with selected model

confirmed.out.train <- glm(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
                       + koi_period + koi_depth + koi_prad + koi_teq + koi_model_snr + koi_steff
                       + koi_kepmag, family = binomial, data=train_set)

trainval.prob <- predict(confirmed.out.train, validation_set, type="response")

# probability chosen with coords(roc.out, "best")
trainval.pred <- (trainval.prob > 0.47) + 0


# Confusion matrix
conf.matrix.trainval<- table(validation_set$confirmed, trainval.pred)
dimnames(conf.matrix.trainval) <- list("Predictions" = c(0, 1), "Confirmed" = c(0, 1))
conf.matrix.trainval

# Performance measures
perf.measure(validation_set$confirmed, trainval.pred)

# ROC curve
roc.out <- roc(validation_set$confirmed, trainval.prob)
plot(roc.out, print.auc=TRUE, legacy.axes=TRUE, xlab="False positive rate", ylab="True positive rate")

# Our chosen model performs well also on the validation set, with an AUC of 0.920
# It is less precise on unseen data, but recall improved.


# LDA

lda.fit <- lda(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
               + koi_fpflag_ec + koi_period + koi_impact + koi_duration
               + koi_depth + koi_prad + koi_teq + koi_insol
               + koi_model_snr + koi_steff + koi_slogg + koi_srad
               + koi_kepmag, data=train_set)
lda.fit

# The graph is not very clear as the density of the observations is very sparse.
plot(lda.fit)

# Predicting on the validation set
lda.pred <- predict(lda.fit, validation_set)

lda.pred$posterior[10:20,] # we check the probabilities for an observation
                           # to be part of each class
lda.pred$class[10:20] # we check the actual predictions

lda.class <- rep(0, 2399) # initializing the vectors of 0s
lda.class[lda.pred$posterior[,2]>= 0.5] <- 1

# For the default threshold (0.5) we should get 0
sum(lda.class!=lda.pred$class)

perf.measure(validation_set$confirmed, lda.pred$class, lab.pos = 1)

# We obtained really good recall, lower precision and specificity than
# the logistic regression.


# QDA

qda.fit <- qda(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
               + koi_fpflag_ec + koi_period + koi_impact + koi_duration
               + koi_depth + koi_prad + koi_teq + koi_insol
               + koi_model_snr + koi_steff + koi_slogg + koi_srad
               + koi_kepmag, data=train_set)
qda.fit
# We get the error "rank deficiency in group 1".

vif_model <- lm(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
                + koi_fpflag_ec + koi_period + koi_impact + koi_duration
                + koi_insol
                + koi_steff +
                + koi_kepmag, data=train_set)
vif(vif_model)

# Increasing the train size (from 0.7 to 0.9)
train_size <- floor(0.9 * nrow(koi.r))

train_indices <- sample(seq_len(nrow(koi.r)), size = train_size)

# Split the data into training and validation sets
train_set <- koi.r[train_indices, ]
validation_set <- koi.r[-train_indices, ]

qda.fit <- qda(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_period 
               + koi_impact + koi_duration + koi_insol+ koi_steff
               + koi_kepmag, data=train_set)
qda.fit

# We were able to solve the error "rank deficiency in group 1" by checking for
# multicollinearity, removing the variables with higher VIF values,
# and incrementing the size of the training set.


qda.class <- predict(qda.fit, validation_set)$class
table(qda.class, validation_set$confirmed)

perf.measure(validation_set$confirmed, qda.class, lab.pos = 1)



# K-Nearest Neighbors (KNN)

set.seed(42)  
n <- nrow(koi.r)
validation_indices <- sample(1:n, size = 0.2 * n)  # 20% validation set
train_indices <- setdiff(1:n, validation_indices)  # remaining 80% training set

train_set <- koi.r[train_indices, ]
validation_set <- koi.r[validation_indices, ]


# Prepare the data for k-NN
train.X <- as.matrix(train_set[, c("koi_fpflag_nt", "koi_fpflag_ss", "koi_fpflag_co", 
                                   "koi_fpflag_ec", "koi_period", "koi_impact", 
                                   "koi_duration", "koi_depth", "koi_prad", 
                                   "koi_teq", "koi_insol", "koi_model_snr", 
                                   "koi_steff", "koi_slogg", "koi_srad", 
                                   "koi_kepmag")])
validation.X <- as.matrix(validation_set[, c("koi_fpflag_nt", "koi_fpflag_ss", "koi_fpflag_co", 
                                             "koi_fpflag_ec", "koi_period", "koi_impact", 
                                             "koi_duration", "koi_depth", "koi_prad", 
                                             "koi_teq", "koi_insol", "koi_model_snr", 
                                             "koi_steff", "koi_slogg", "koi_srad", 
                                             "koi_kepmag")])
train.confirmed <- train_set$confirmed
validation.confirmed <- validation_set$confirmed

# k-NN with k = 1

knn.pred <- knn(train.X, validation.X, train.confirmed, k = 1)
conf.matrix <- table(knn.pred, validation.confirmed)
conf.matrix

perf.measure(validation.confirmed, knn.pred, lab.pos = 1)

# k-NN with k = 3

knn.pred <- knn(train.X, validation.X, train.confirmed, k = 3)
conf.matrix <- table(knn.pred, validation.confirmed)
conf.matrix

perf.measure(validation.confirmed, knn.pred, lab.pos = 1)

# k-NN with k = 11

knn.pred <- knn(train.X, validation.X, train.confirmed, k = 11)
conf.matrix <- table(knn.pred, validation.confirmed)
conf.matrix

perf.measure(validation.confirmed, knn.pred, lab.pos = 1)

# k-NN with k = 51

knn.pred <- knn(train.X, validation.X, train.confirmed, k = 51)
conf.matrix <- table(knn.pred, validation.confirmed)
conf.matrix

perf.measure(validation.confirmed, knn.pred, lab.pos = 1)

# k-NN with k = 71

knn.pred <- knn(train.X, validation.X, train.confirmed, k = 71)
conf.matrix <- table(knn.pred, validation.confirmed)
conf.matrix

perf.measure(validation.confirmed, knn.pred, lab.pos = 1)

# We saw that K higher than 71 did not improve the overall error rate significantly,
# and we stopped testing for higher values. We choose odd values for k to avoid ties,
# as we have a binary classification problem.
# K-NN is the worst method tested so far. Lots of data, combined with many predictors
# cause underfitting.


# LOOCV

glm.fit <- glm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                 koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                 koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag, 
               data = koi.r, family = gaussian)

# cv.err <- cv.glm(koi.r, glm.fit)

# Number of groups (should be equal to the number of observations for LOOCV)
print(cv.err$K)

# Cross-validation error (first component: raw estimate, second: adjusted)
cv.err$delta

# LOOCV is too computationally expensive for our dataset.




regfit.full <- regsubsets(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                            koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                            koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag, data=koi.r)


reg.summary <- summary(regfit.full)
reg.summary$outmat

names(reg.summary)
reg.summary$rsq
reg.summary$bic


# residual sum of squares
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="l")

# adjusted-R^2 with its largest value
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted Rsq",type="l")
i <- which.max(reg.summary$adjr2)
points(i,reg.summary$adjr2[i], col="red",cex=2,pch=20)
text(i,reg.summary$adjr2[i], i, pos=1)

# Mallow's Cp with its smallest value
plot(reg.summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
i <- which.min(reg.summary$cp)
points(i,reg.summary$cp[i],col="red",cex=2,pch=20)
text(i,reg.summary$cp[i], i, pos=3)

# BIC with its smallest value
plot(reg.summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
i <- which.min(reg.summary$bic)
points(i,reg.summary$bic[i],col="red",cex=2,pch=20)
text(i,reg.summary$bic[i], i, pos=3)

# We need to return models with more variables
regfit.full <- regsubsets(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                            koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                            koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                            data=koi.r, nvmax=16)

reg.summary <- summary(regfit.full)
reg.summary$outmat

names(reg.summary)
reg.summary$rsq
reg.summary$bic


# residual sum of squares
plot(reg.summary$rss,xlab="Number of Variables",ylab="RSS",type="l")

# adjusted-R^2 with its largest value
plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted Rsq",type="l")
i <- which.max(reg.summary$adjr2)
points(i,reg.summary$adjr2[i], col="red",cex=2,pch=20)
text(i,reg.summary$adjr2[i], i, pos=1)

# Mallow's Cp with its smallest value
plot(reg.summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
i <- which.min(reg.summary$cp)
points(i,reg.summary$cp[i],col="red",cex=2,pch=20)
text(i,reg.summary$cp[i], i, pos=3)

# BIC with its smallest value
plot(reg.summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
i <- which.min(reg.summary$bic)
points(i,reg.summary$bic[i],col="red",cex=2,pch=20)
text(i,reg.summary$bic[i], i, pos=3)


# The best model, according to BIC, is the one with 13 predictors.

coef(regfit.full, 13)

best.bic <- lm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co 
               + koi_fpflag_ec + koi_period + koi_impact + koi_duration
               + koi_depth + koi_prad + koi_teq + koi_insol + koi_model_snr
               + koi_slogg, data=koi.r)
summary(best.bic)

# Check how R-squared changes 

reg.summary$rsq

# Diagnostic plots

plot(best.bic, which = 1)
plot(best.bic, which = 2)
plot(best.bic, which = 3)
plot(best.bic, which = 5)


transformed.fit <- lm(log(koi_score + 1e-8 / (1 - koi_score + 1e-8)) ~ 
                        koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                        log(koi_period + 1) + log(koi_impact + 1) + log(koi_duration + 1) + 
                        log(koi_depth + 1) + log(koi_prad + 1) + log(koi_teq + 1) + 
                        log(koi_insol + 1) + log(koi_model_snr + 1) + log(koi_slogg + 1), 
                      data = koi.r)

plot(transformed.fit, which = 1)
plot(transformed.fit, which = 2)
plot(transformed.fit, which = 3)
plot(transformed.fit, which = 5)


transformed.fit <- lm(log(koi_score + 1e-8 / (1 - koi_score + 1e-8)) ~ 
                        koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                        poly(koi_period + 1, 2) + poly(koi_duration + 1, 2) + 
                        poly(koi_depth + 1, 2) + poly(koi_prad + 1, 2) + poly(koi_teq + 1, 2) + 
                        poly(koi_insol + 1, 2) + poly(koi_impact + 1, 2) + poly(koi_model_snr + 1, 2) +
                        koi_period:koi_depth + koi_prad:koi_teq + koi_period:koi_insol + 
                        koi_duration:koi_impact + koi_model_snr:koi_insol, 
                      data = koi.r)


transformed.fit <- regsubsets(koi_score ~ 
                        log(koi_fpflag_nt + 1e-8) + log(koi_fpflag_ss + 1e-8) + log(koi_fpflag_co + 1e-8) + log(koi_fpflag_ec + 1e-8) + 
                        poly(koi_period, 10) + log(koi_period + 1) +
                        poly(koi_duration, 10) + log(koi_duration + 1) +
                        poly(koi_depth, 10) + log(koi_depth + 1) +
                        poly(koi_prad, 5) + log(koi_prad + 1) +
                        poly(koi_teq, 10) + log(koi_teq + 1) +
                        poly(koi_insol, 10) + log(koi_insol + 1) +
                        poly(koi_impact, 10) + log(koi_impact + 1) +
                        poly(koi_model_snr, 10) + log(koi_model_snr + 1) +
                        koi_period:koi_depth + koi_prad:koi_teq + koi_period:koi_insol + 
                        koi_duration:koi_impact + koi_model_snr:koi_insol, 
                      data = koi.r, nvmax = 25)
# Error: Exhaustive search will be S L O W, must specify really.big=T

# We tried to solve non-linearity and heteroschedasticity issues, but nothing
# seems to work. With 13 predictors there are too many combination to test.
# After a few tries we stopped. The smartest way to fit better this dataset wuold
# be to use regressions that are robust to heteroschedasticity, and/or non-linear
# regression methods.

# With regsubsets we selected the model with 13 predictors based on BIC.
# Let's try other methods.

#Forward Inclusion

regfit.fwd <- regsubsets(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                           koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                           koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                         data=koi.r, nvmax=16, method="forward")
plot(regfit.fwd, scale="bic")

fwd.summary <- summary(regfit.fwd)
plot(fwd.summary$bic, xlab = "Number of Variables", ylab="BIC",type='l')
i <- which.min(fwd.summary$bic)
i
points(i,reg.summary$bic[i],col="red",cex=2,pch=20)

# We get, as before to 13 predictors, based on BIC.

# Backward elimination

regfit.bwd <- regsubsets(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                           koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                           koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                         data=koi.r, nvmax=16, method="backward")
plot(regfit.bwd, scale="bic")

bwd.summary <- summary(regfit.bwd)
plot(bwd.summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
i <- which.min(bwd.summary$bic)
i
points(i,bwd.summary$bic[i],col="red",cex=2,pch=20)

# Same with the backward approach

coef(regfit.full, 13)
coef(regfit.fwd, 13)
coef(regfit.bwd, 13)

# Coefficients are the same for every approach.


# Stepwise using the step function

n <- nrow(koi.r)

mod.F <- lm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
              koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
              koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag, data=koi.r)
mod.bwd.bic <- step(mod.F, direction="backward", k=log(n), trace=1, steps=1000)
summary(mod.bwd.bic)

# Same coefficients as before


# K-fold Cross Validation

# Set a value for k (number of folds)
k <- 10

# Associate a fold number to each observation
set.seed(42)
# set.seed(222)
# set.seed(123)
# set.seed(231)

folds <- sample(1:k, nrow(koi.r), replace = TRUE)
folds[1:20]
table(folds)

# For each of the k=10 folds we are going to compute p=16 validation errors
# (since we have 16 predictors) and thus we need a 10x16 matrix
cv.errors <- matrix(NA, k, 16)
colnames(cv.errors) <- 1:16

# Apply the validation set approach k times, one for each fold
for (j in 1:k) {
  best.fit <- regsubsets(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                           koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                           koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                         data = koi.r[folds != j,], nvmax = 16)
  test.mat <- model.matrix(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                             koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                             koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                           data = koi.r[folds == j,])
  for (i in 1:16) {
    beta <- coef(best.fit, id = i)
    pred <- test.mat[, names(beta)] %*% beta
    cv.errors[j, i] <- mean((koi.r$koi_score[folds == j] - pred)^2)
  }
}

# Compute the mean of CV errors across the k folds
mean.cv.errors <- apply(cv.errors, 2, mean)
mean.cv.errors

# Plot the mean CV errors
plot(mean.cv.errors, type = 'b', xlab = 'Number of Predictors', ylab = 'CV Error')
i <- which.min(mean.cv.errors)
points(i, mean.cv.errors[i], col = "red", cex = 2, pch = 20)

# Fit the selected model on the full dataset
reg.best <- regsubsets(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                         koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                         koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                       data = koi.r, nvmax = 16)
coef(reg.best, i)

# With the K-fold Cross Validation method, we found 11 predictors for the best model.







# Ridge and Lasso Regressions on "confirmed" binary variable

confirmed.out <- glm(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co
                     + koi_fpflag_ec + koi_period + koi_impact + koi_duration
                     + koi_depth + koi_prad + koi_teq + koi_insol
                     + koi_model_snr + koi_steff + koi_slogg + koi_srad
                     + koi_kepmag, family = binomial, data=koi.r)
summary(confirmed.out)



# Prepare the data
# Create the design matrix and response vector
X <- model.matrix(confirmed ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec + 
                    koi_period + koi_impact + koi_duration + koi_depth + koi_prad + koi_teq + 
                    koi_insol + koi_model_snr + koi_steff + koi_slogg + koi_srad + koi_kepmag,
                  data = koi.r)[,-1]
y <- koi.r$confirmed

# Split data into training and test sets
set.seed(42)
train <- sample(1:nrow(X), nrow(X) / 2)
test <- setdiff(1:nrow(X), train)

# Define the grid of lambda values
grid <- 10^seq(10, -2, length = 100)

###########################
# RIDGE REGRESSION
###########################

# Fit ridge regression model on training set
ridge.mod <- glmnet(X[train, ], y[train], alpha = 0, lambda = grid, family = "binomial")

# Plot the ridge regression coefficients
plot(ridge.mod, xvar = "lambda", label = TRUE)

# Perform cross-validation to choose the best lambda
set.seed(42)
cv.out <- cv.glmnet(X[train, ], y[train], alpha = 0, lambda = grid, nfolds = 10, family = "binomial")
plot(cv.out)

# Identify the best lambda
bestlam.ridge <- cv.out$lambda.min
bestlam.ridge

# Predict on the test set using the best lambda
ridge.pred <- predict(ridge.mod, s = bestlam.ridge, newx = X[test, ], type = "response")

# Calculate the test MSE
ridge.mse <- mean((ridge.pred - y[test])^2)
ridge.mse

# Fit the final model on the full dataset
ridge.final <- glmnet(X, y, alpha = 0, lambda = grid, family = "binomial")
coef(ridge.final, s = bestlam.ridge)

###########################
# LASSO REGRESSION
###########################

# Fit lasso regression model on training set
lasso.mod <- glmnet(X[train, ], y[train], alpha = 1, lambda = grid, family = "binomial")

# Plot the lasso regression coefficients
plot(lasso.mod, xvar = "lambda", label = TRUE)

# Perform cross-validation to choose the best lambda
set.seed(42)
cv.out.lasso <- cv.glmnet(X[train, ], y[train], alpha = 1, lambda = grid, nfolds = 10, family = "binomial")
plot(cv.out.lasso)

# Identify the best lambda
bestlam.lasso <- cv.out.lasso$lambda.min
bestlam.lasso

# Predict on the test set using the best lambda
lasso.pred <- predict(lasso.mod, s = bestlam.lasso, newx = X[test, ], type = "response")

# Calculate the test MSE
lasso.mse <- mean((lasso.pred - y[test])^2)
lasso.mse

# Fit the final model on the full dataset
lasso.final <- glmnet(X, y, alpha = 1, lambda = grid, family = "binomial")
coef(lasso.final, s = bestlam.lasso)

# Compare test MSE of ridge and lasso
ridge.mse
lasso.mse

# Print coefficients of the final lasso model
print(coef(lasso.final, s = bestlam.lasso))






































############## OUTDATED
# Regression only on the categorical variables
mod.out <- glm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_fpflag_ec, family = binomial, data = koi.r)
summary(mod.out)

mod.out.resid <- residuals(mod.out, type="deviance")
plot(mod.out.resid~fitted(mod.out))
qqnorm(mod.out.resid)
qqline(mod.out.resid)


# Regression only on the categorical variables - 1
mod.outR <- glm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co, family = binomial, data = koi.r)
summary(mod.outR)

mod.outR.resid <- residuals(mod.outR, type="deviance")
plot(mod.outR.resid~fitted(mod.outR))
qqnorm(mod.outR.resid)
qqline(mod.outR.resid)


# One dummy's information is completely contained in the combination
# of the other 3 dummies?

# The full model has a lower residual variance, and lower AIC.
# At the same time, koi_fpflag_ec has a p-value of 0.381,
# which tells us that the coefficient is not statistically significant
# thus indicating that koi_fpflag_ec does not have a strong predictive effect
# on koi_score in this model.

# Should we use the full model or the reduced one?


# anova test
anova(mod.outR, mod.out, test="Chisq")

# The very low p-value (< 2.2e-16) indicates that the addition of koi_fpflag_ec
# to the model significantly improves the model fit.
# The reduction in residual deviance from 3119.7 (Model 1) to 2591.5 (Model 2)
# also shows that Model 2 fits the data better than Model 1.

# Even though the coefficient for koi_fpflag_ec in the individual model output
# was not significant, when comparing models, the addition of this predictor
# does significantly improve the model fit overall.
# This is probably due to interdependencies between predictors, and we know
# for sure these interdependencies exists since all 4 variables are dummies,
# and this means that one of them is completely explained by the other three.





## Regression on all variables
mod.out <- glm(koi_score ~ ., family = binomial, data = koi.r)
summary(mod.out)
# Residual deviance: 2240.3, AIC: 2918.3

# How can we improve our fit?

##  Feature selection
# Let's start by simply removing all non-significant predictors.

mod.out.reduced <- glm(koi_score ~ koi_fpflag_nt + koi_fpflag_ss + koi_fpflag_co + koi_period + 
                         koi_impact + koi_duration + koi_depth + koi_teq + koi_insol, 
                       family = binomial, data = koi.r)
summary(mod.out.reduced)
# Residual deviance: 2782.9, AIC: 3475.9

# anova test
anova(mod.out.reduced, mod.out, test="Chisq")
# Pr(>Chi): 2.2e-16 ***

# We can clearly see how our fit got worse, and the anova test confirms it.


## Including interaction terms
mod.out.interactions <- glm(koi_score ~ koi_fpflag_nt * koi_fpflag_ss + koi_fpflag_nt * koi_fpflag_co + 
                              koi_fpflag_ss * koi_fpflag_co + koi_period + koi_impact + 
                              koi_duration + koi_depth + koi_teq + koi_insol, 
                            family = binomial, data = koi.r)
summary(mod.out.interactions)


## Including polynomial terms
mod.out.polynomial <- glm(koi_score ~ poly(koi_period, 2) + poly(koi_impact, 2) + 
                            poly(koi_duration, 2) + poly(koi_depth, 2) + koi_fpflag_nt + 
                            koi_fpflag_ss + koi_fpflag_co + koi_teq + koi_insol, 
                          family = binomial, data = koi.r)
summary(mod.out.polynomial)


## Fitting the model with transformed features
mod.out.transformed <- glm(koi_score ~ log(koi_period) + sqrt(koi_duration) + koi_fpflag_nt + 
                             koi_fpflag_ss + koi_fpflag_co + koi_impact + koi_depth + 
                             koi_teq + koi_insol, family = binomial, data = koi.r)
summary(mod.out.transformed)



