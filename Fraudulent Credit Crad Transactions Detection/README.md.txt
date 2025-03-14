# 🚨 Fraudulent Credit Card Transactions Detection

## 📌 Overview
This project focuses on **real-time credit card fraud detection** using **deep learning**. The challenge lies in the **highly imbalanced dataset**, where fraudulent transactions are extremely rare. By leveraging **SMOTE**, feature engineering, and a deep learning model, we aim to improve fraud detection accuracy.

## 📊 Dataset
- **Source**: Kaggle - European Credit Card Transactions Dataset
- **Features**:
  - **28 anonymized transaction features (V1 to V28)**
  - **Time & Amount** (normalized)
  - **Class**: **0 (legitimate)**, **1 (fraudulent)**
- **Imbalance**: Only **0.17% of transactions** are fraudulent.

## 🔍 Key Analyses
1. **Exploratory Data Analysis (EDA)**:
   - **Class Imbalance Visualization**: Highlighted fraud distribution.
   - **Transaction Time & Amount Analysis**: Identified trends in fraudulent transactions.
   - **Correlation Heatmap**: Detected key influencing variables.

2. **Feature Engineering**:
   - Normalized **Time & Amount**.
   - Applied **PCA** for dimensionality reduction.

3. **Handling Class Imbalance**:
   - Used **SMOTE (Synthetic Minority Over-sampling Technique)** to generate synthetic fraud cases.

4. **Deep Learning Model**:
   - **Architecture**:
     - **Input Layer**: 64 neurons
     - **Hidden Layers**: Two layers with **32 neurons each**, ReLU activation.
     - **Dropout (20%)**: Prevents overfitting.
     - **Output Layer**: Sigmoid activation for binary classification.

5. **Evaluation Metrics**:
   - **Precision-Recall Curve (AUPRC)**: Achieved **0.83**, indicating strong fraud detection capability.
   - **Confusion Matrix**: High **true positive rate**, minimal false alarms.
   - **F1-Score**: Balanced fraud detection performance.

## 📈 Visualizations
- **Class Distribution Heatmap**: Shows fraud cases vs. legitimate transactions.
- **Feature Importance Plot**: Highlights key predictors of fraudulent activity.
- **Confusion Matrix**: Visualizes model performance.

## 📌 Key Findings
- **Fraudulent transactions tend to be of lower amounts** to avoid detection.
- **Transaction time correlates with fraud risk**—certain time windows show higher fraud activity.
- **Feature engineering significantly improves model performance**.

## 🛠 Tools Used
- **Python** (Pandas, NumPy, Seaborn, Matplotlib)
- **Scikit-learn** (SMOTE, PCA, Feature Engineering)
- **TensorFlow/Keras** (Deep Learning Model)
- **Matplotlib & Seaborn** (Visualizations)

## 🎯 Future Work
- Experiment with **Random Forest & XGBoost** for comparison.
- Implement **real-time fraud detection pipeline**.
- Further optimize **hyperparameters & dropout rates**.

---
