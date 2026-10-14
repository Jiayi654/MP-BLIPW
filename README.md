# MP-BLIPW

This repository contains the R code accompanying the manuscript on Multiphase Best Linear Inverse Probability Weighting (MP-BLIPW) for regression analysis with missing covariates under multiphase sampling designs.

## Repository Structure

### Simulation Study

The simulation code evaluates the performance of the following estimators:

- Single Inverse Probability Weighting (SIPW)
- Best Linear Inverse Probability Weighting (BLIPW)
- Multiphase Best Linear Inverse Probability Weighting (MP-BLIPW)

The simulation studies compare estimation bias, efficiency, standard errors, and confidence interval coverage under different multiphase sampling scenarios.

### NHANES Application

The empirical analysis applies MP-BLIPW to NHANES 2017–2018 data to investigate the association between high-density lipoprotein (HDL) cholesterol and DEXA-measured total body fat percentage.

Body mass index (BMI) and hip circumference are incorporated as auxiliary variables in the multiphase sampling design.

### Software

All simulations and data analyses were conducted using R.

### Reproducibility

The repository contains R scripts for implementing the proposed estimator, conducting simulation studies, and performing the NHANES application.

Run the R scripts in numerical order, ensuring that all required packages and data are available.

## Authors

- Osvaldo Espin Garcia
- Jiayi Wang
- Kathleen Zang
