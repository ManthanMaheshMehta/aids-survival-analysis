# Survival Prediction in AIDS Patients Using Traditional and Machine Learning-Based Survival Models

This project compares traditional survival models with machine learning-based approaches to predict survival outcomes among AIDS patients. It uses data from the ACTG Study 175, a clinical trial investigating the effectiveness of antiretroviral therapies.

## 📊 Objective
To evaluate and compare the predictive performance of:
- **Cox Proportional Hazards Model**
- **Random Survival Forest (RSF)**
- **Gradient Boosting Machine (GBM)**

## 📁 Dataset
- **Source**: ACTG Study 175 (Hammer et al., 1996)
- **Sample Size**: 2,139 patients
- **Outcome**: Time to event (death or censoring)

## 🧪 Methods
- Time-to-event analysis using `survival`, `randomForestSRC`, and `gbm` packages in R
- Risk score calculation based on Cox model
- Evaluation metrics: **Concordance Index (C-index)**, **CRPS**, **Prediction Error**
- Variable importance ranking
- Visualization: Kaplan-Meier curves, Manhattan plots, risk score plots

## 🚀 Tools Used
- R, R Markdown
- Packages: `survival`, `survminer`, `randomForestSRC`, `gbm`, `ggplot2`
- Git & GitHub
- R Shiny (for dashboard prototype)

## 📈 Key Results
- RSF showed slightly better performance (C-index = 0.616) compared to Cox (0.614)
- CD4 counts, Karnofsky score, and treatment history were top predictors
- GBM underperformed due to overfitting in this dataset

## 📂 Files
- `data/` – cleaned dataset  
- `code/` – analysis scripts (Cox, RSF, GBM, visualizations)  
- `figures/` – plots and charts used in the report  
- `report/` – final write-up and presentation slides  
- `shiny_app/` – interactive dashboard (optional)

## 📚 Reference
Hammer et al., 1996. A randomized, controlled trial of four antiretroviral treatments in HIV-infected patients.

## 📬 Contact
Manthan Mahesh Mehta  
📧 mehtmanthan1999@gmail.com 

