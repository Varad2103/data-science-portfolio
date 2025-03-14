# 💳 CreditWise: A Comprehensive Analysis of Credit Scores

## 📌 Overview
This project explores **credit score analysis** using statistical and probabilistic methods. It applies **Bayesian theorem, joint probability distributions, and factor analysis** to uncover key determinants affecting an individual's credit score. The study provides actionable insights for financial institutions and individuals to improve creditworthiness.

## 📊 Dataset
- **Source**: [Kaggle - Credit Score Classification Dataset](https://www.kaggle.com/datasets/sujithmandala/credit-score-classification-dataset/data)
- **Features**:
  - **Demographics**: Age, Gender, Marital Status, Number of Children
  - **Financial Factors**: Income, Home Ownership, Education
  - **Target Variable**: Credit Score (Low, Average, High)
- **Preprocessing**: No missing values, cleaned dataset ready for analysis.

## 🔍 Key Analyses
1. **Bayesian Analysis**:
   - Evaluated conditional probabilities of credit scores based on demographic attributes.
   - Found that **homeownership** and **income levels** significantly impact credit scores.
  
2. **Joint Probability Distributions**:
   - Measured relationships between **age, income, and education** on creditworthiness.
   - **High income & homeownership** lead to better credit scores.

3. **Factor Analysis**:
   - Used **KMO test & eigenvalues** to identify latent factors.
   - **Age and Income** strongly correlated, influencing overall creditworthiness.

## 📈 Visualizations
- **Correlation Heatmap**: Identifies key attributes affecting credit scores.
- **Probability Distributions**: Displays likelihood of different credit scores.
- **Outlier Detection**: Used **IQR & Z-score methods** to detect anomalies.
- **Factor Loading Scree Plot**: Shows variance explained by key financial indicators.

## 📌 Key Findings
- **Females are more likely to have lower credit scores** compared to males.
- **Married individuals tend to have higher credit scores**.
- **Homeownership and stable income** are strong indicators of high creditworthiness.
- **Young individuals (<30 years) have lower probabilities of high credit scores**.
- **Education plays a crucial role**—higher education levels correlate with better credit scores.

## 🛠 Tools Used
- **Python** (NumPy, Pandas, Seaborn, Matplotlib)
- **Scipy** (Statistical Testing: Bayes’ Theorem, ANOVA)
- **FactorAnalyzer** (Factor Analysis)
- **Sklearn** (Outlier detection, Standardization)

## 🎯 Recommendations
- Banks should prioritize **income, homeownership, and education** in credit score assessments.
- Young individuals should focus on **building financial stability** before applying for loans.
- Financial institutions can **adjust risk models** based on demographic patterns.

---
