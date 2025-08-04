# NASA Exoplanet Dataset Analysis: Predicting KOI Probabilities and Classification

## Project Overview

This project analyzes the NASA Exoplanet Dataset to predict probabilities and classify Kepler Objects of Interest (KOIs). KOIs are stars suspected of hosting one or more transiting planets, detected by NASA's Kepler mission. The analysis focuses on predicting the `koi_score` (likelihood of being a true exoplanet) and classifying KOIs as confirmed planets or false positives.

## Dataset Description

- **Source**: NASA Exoplanet Archive (KOI Table - Cumulative List)
- **Size**: 9,564 observations with 27 variables
- **Format**: CSV file
- **Each observation**: Represents a different KOI event
- **Unique stars**: 8,214 stars (some have multiple KOI events)

### Key Variables

| Variable | Description |
|----------|-------------|
| `koi_score` | Quantitative measure of likelihood that a KOI is a true exoplanet (0-1) |
| `koi_disposition` | Final classification: CONFIRMED, CANDIDATE, or FALSE POSITIVE |
| `koi_fpflag_*` | Binary flags indicating different types of false positive signals |
| `koi_period` | Orbital period of the potential exoplanet (days) |
| `koi_depth` | Transit depth (amount of brightness decrease) |
| `koi_duration` | Duration of transit event |
| `koi_prad` | Estimated planet radius (Earth radii) |
| `koi_teq` | Estimated equilibrium temperature (Kelvin) |
| Stellar properties | Host star characteristics (temperature, radius, magnitude, etc.) |

## Project Structure

```
├── koi.R                # Main analysis script
├── report.Rmd           # R Markdown report
├── koi.csv              # Original dataset
└── README.md            
```

## Analysis Pipeline

### 1. Data Cleaning and Preprocessing
- **Duplicate handling**: No duplicate KOIs found, but 1,350 stars have multiple KOI events
- **Missing values**: Systematic handling of NAs and ghost values
- **Data validation**: Fixed erroneous values (e.g., 465 → 0 in binary flag)
- **Variable transformation**: Log transforms, unit conversions (Kelvin to Celsius, ppm to percentages)
- **Feature selection**: Removed non-predictive variables (coordinates, identifiers, timestamps)

### 2. Exploratory Data Analysis (EDA)
- Distribution analysis of all variables
- Correlation matrix visualization
- Relationship exploration between predictors and response variables
- Identification of outliers and data patterns

### 3. Multiple Linear Regression
**Objective**: Predict `koi_score` using various KOI characteristics

**Key Findings**:
- R² ≈ 0.75 with all predictors
- Four binary false positive flags (`koi_fpflag_*`) explain ~67% of variance
- Continuous variables show limited linear relationship with `koi_score`
- Diagnostic plots reveal heteroscedasticity and non-linear relationships
- Transformations (log, linear-log, log-log) didn't significantly improve model fit

### 4. Model Selection Techniques
- **Best Subset Selection**: 13 predictors optimal (BIC criterion)
- **Forward Selection**: Same 13-predictor model
- **Backward Elimination**: Consistent results
- **K-Fold Cross Validation**: 11 predictors optimal

### 5. Classification Analysis
**Objective**: Classify KOIs as confirmed planets vs. not confirmed

#### Logistic Regression
- **Primary model**: Predicting `confirmed` status
- **Performance**: AUC = 0.920 on validation set
- **Feature selection**: Constraint-based backward elimination
- **Optimal threshold**: 0.44 (maximizing sensitivity + specificity)

#### Alternative Classification Methods
- **LDA**: Poor performance due to violated normality assumptions
- **QDA**: High recall but poor precision and specificity  
- **K-NN**: Underperformed across all k values tested
- **Ridge & Lasso Regression**: Equivalent MSE performance, Lasso preferred for sparsity

## Key Results

### Multiple Linear Regression on `koi_score`
- **Best model R²**: ~0.75
- **Primary predictors**: Binary false positive flags
- **Challenge**: Non-linear relationships and heteroscedasticity

### Classification Performance (Logistic Regression)
- **AUC**: 0.920
- **Precision**: Moderate (varies by threshold)
- **Recall**: High performance
- **Best approach**: Regularized logistic regression

## Dependencies

```r
library(car)        # VIF calculation
library(pROC)       # ROC curves
library(MASS)       # LDA/QDA
library(class)      # K-NN
library(boot)       # Cross-validation
library(leaps)      # Best subset selection
library(glmnet)     # Ridge/Lasso regression
library(ISLR2)      # Datasets and functions
library(corrplot)   # Correlation visualization
```

## Usage

1. **Load the data**:
   ```r
   koi <- read.csv("koi.csv")
   ```

2. **Run the complete analysis**:
   ```r
   source("koi.R")
   ```

3. **Generate the report**:
   ```r
   rmarkdown::render("report.Rmd")
   ```

## Key Insights

1. **False positive flags are the strongest predictors** of both `koi_score` and confirmation status
2. **Continuous variables** (period, depth, duration, etc.) provide information about classification certainty but not classification direction
3. **Non-linear relationships** dominate the dataset, limiting linear model effectiveness
4. **Regularization techniques** show promise for handling the high-dimensional feature space
5. **Classification is more successful** than regression prediction of exact scores

## Future Improvements

- **Non-linear models**: Tree-based methods, neural networks, or GAMs
- **Feature engineering**: Domain-specific transformations and interactions
- **Ensemble methods**: Combining multiple model predictions
- **Deep learning**: For capturing complex non-linear patterns
- **Domain expertise**: Incorporating astrophysics knowledge for better feature selection

## Authors

- Filippo Forcella
- Angelo Guarino

## License

This project uses publicly available NASA data. Please cite the NASA Exoplanet Archive when using this analysis.

## References

- NASA Exoplanet Archive: https://exoplanetarchive.ipac.caltech.edu/
- Kepler Mission: https://www.nasa.gov/mission_pages/kepler/main/index.html