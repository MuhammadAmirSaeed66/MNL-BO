# 🧠 MNL-BO: Multinomial Logit Bayesian Optimization

<p align="center">
  <img src="https://github.com/MuhammadAmirSaeed66.png" width="120" style="border-radius:50%">
</p>

<p align="center">
  <b>Reproducible R Implementation for Bayesian Optimization with Categorical and Mixed Variables Using a Multinomial Logit Surrogate</b>
</p>

<p align="center">
  <a href="https://github.com/MuhammadAmirSaeed66/MNL-BO/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg">
  </a>
  <!-- TODO: Replace this Zenodo badge after creating a Zenodo DOI -->
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

Bayesian Optimization (BO) is widely used for optimizing expensive black-box functions. Traditional BO methods rely on Gaussian process surrogates, which perform well in continuous domains but encounter difficulties in categorical and mixed-variable settings.

This repository implements **MNL-BO**, a preference-based Bayesian optimization framework using a **Multinomial Logit (MNL)** surrogate trained from pairwise comparisons.

---

## ✨ Key Contributions

* Multinomial Logit surrogate for BO
* Natural handling of categorical variables
* Unified framework for mixed-variable optimization
* Preference-based learning using pairwise comparisons
* Evaluated on:

  * Categorical benchmark
  * Traveling Salesman Problem
  * Pressure vessel design

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
```

---

## ⚙️ Installation

```bash
git clone https://github.com/MuhammadAmirSaeed66/MNL-BO.git
cd MNL-BO
```

```r
install.packages(c("tidyverse","data.table","mlogit","ggplot2"))
```

---

## ▶️ Usage

Run all experiments:

```r
source("scripts/run_experiments.R")
```

Run individual studies:

```r
source("scripts/categorical_case.R")
source("scripts/tsp_case.R")
source("scripts/pressure_vessel_case.R")
```

---

## 🧪 Reproducibility

All results in the manuscript can be reproduced using:

```r
source("scripts/run_experiments.R")
```

---

## 📚 Citation

If you use this code, please cite:

```bibtex
@article{saeed2026mnlbo,
  title   = {Bayesian Optimization for Categorical and Mixed Variables Using a Multinomial Logit Surrogate},
  author  = {Saeed, Muhammad Amir and others},
  journal = {Algorithms},
  year    = {2026}
}
```

---

## 👨‍💻 Author

**Muhammad Amir Saeed**

<p>
  <a href="https://github.com/MuhammadAmirSaeed66">
    <img src="https://img.shields.io/badge/GitHub-MuhammadAmirSaeed66-black?logo=github">
  </a>
</p>

---

## 📜 License

This repository is licensed under the MIT License.

---

## ⭐ Support

If you find this repository useful, please consider giving it a star ⭐

---
