# 🧠 MNL-BO: Multinomial Logit Bayesian Optimization

<p align="center">
  <img src="https://github.com/MuhammadAmirSaeed66.png" width="120" style="border-radius:50%">
</p>

<p align="center">
  <b>Reproducible R Implementation for Bayesian Optimization with Categorical and Mixed Variables</b>
</p>

<p align="center">
  <a href="https://github.com/MuhammadAmirSaeed66/MNL-BO/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg">
  </a>
  <!-- TODO: Replace this Zenodo badge after creating a Zenodo DOI -->
  <!-- <a href="https://zenodo.org/badge/latestdoi/XXXX">
    <img src="https://zenodo.org/badge/XXXX.svg" alt="DOI">
  </a> -->
</p>

---

## 📄 Paper

**Title:**  
*Bayesian Optimization for Categorical and Mixed Variables Using a Multinomial Logit Surrogate*

**Journal:**  
*Algorithms* — MDPI

This repository contains the complete **R source code**, supporting data, scripts, and figures required to reproduce the computational experiments reported in the manuscript.

---

## 🚀 Overview

Bayesian Optimization (BO) is a widely used framework for optimizing expensive black-box functions. Traditional BO methods commonly rely on Gaussian process surrogates, which work well for continuous domains but may struggle with categorical and mixed-variable search spaces.

This repository implements **MNL-BO**, a preference-based Bayesian optimization framework that replaces the Gaussian process surrogate with a **Multinomial Logit (MNL)** model trained from pairwise preference comparisons.

The proposed framework provides an interpretable surrogate model for categorical alternatives while allowing continuous, discrete, and categorical variables to be handled within a unified optimization workflow.

---

## ✨ Key Contributions

- Preference-based Bayesian optimization using a Multinomial Logit surrogate
- Natural handling of categorical and mixed-variable decision spaces
- Pairwise comparison-based learning mechanism
- Acquisition functions balancing exploration and exploitation
- Reproducible experiments for:
  - Purely categorical benchmark problem
  - Combinatorial Traveling Salesman Problem
  - Mixed-variable pressure vessel design optimization
- Comparisons against random search, local search, classical metaheuristics, and SMAC-inspired tree-based BO baselines

---

## 📂 Repository Structure

```bash
MNL-BO/
│── data/              # Supporting datasets
│── R/                 # Core R functions
│── scripts/           # Experiment scripts
│── results/           # Numerical outputs
│── figures/           # Generated figures and demo GIF
│── README.md
│── LICENSE
│── .gitignore
