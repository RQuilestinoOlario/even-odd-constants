# even-odd-constants
Parity analysis of 10 million decimal digits in mathematical constants (π, e, φ, √2, √3, ln 2, ln 10, Catalan G, Euler γ, and lemniscate ϖ).

---

## Overview
This repository contains R code, summary tables, and figures for a large-scale parity and randomness analysis of mathematical constants.  
Each constant’s first 10 million verified digits were obtained from **Alexander J. Yee’s [NumberWorld](https://www.numberworld.org/digits/)** repository.  
Digits were converted to odd/even (1/0) sequences, then tested for balance, local dependence, and randomness using frequency, runs, serial, and approximate-entropy statistics.

---

## Contents
| File | Description |
|------|--------------|
| **01_parity_analysis_10M.R** | Main R script for the full workflow (reading digits, computing parity, tests, and plotting). |
| **summary_10M.csv** | Overall results for each constant: parity proportion, confidence intervals, χ² tests, etc. |
| **kblock_tests_10M.csv** | Serial block-test results for k = 2–5. |
| **digit_counts_10M.csv** | Raw digit frequency counts (0–9) for each constant. |
| **fig01.png** | Cumulative deviation D(n) of odd − even digits. |
| **fig02.png** | Sliding-window parity proportion (w = 10 000 digits). |
| **fig03.png** | Digit-frequency deviations from uniform expectation. |
| **LICENSE** | MIT License. |
| **README.md** | Project overview and usage notes. |

---

## Usage
1. Clone or download this repository.  
2. Obtain digit text files from [NumberWorld](https://www.numberworld.org/digits/) (e.g. `pi.txt`, `e.txt`, …).
3. Edit file paths inside `01_parity_analysis_10M.R` if necessary.  
4. Run in **R (≥ 4.0)** with the following packages installed:
   ```r
   install.packages(c("data.table", "ggplot2"))

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17477614.svg)](https://doi.org/10.5281/zenodo.17477614)
