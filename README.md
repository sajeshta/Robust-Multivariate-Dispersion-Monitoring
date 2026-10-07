# Reproducibility Materials

This repository contains the R code used to reproduce the analyses, numerical results, control limits, Average Run Length (ARL) results, and figures reported in the manuscript:

**“Monitoring Process Variability with Robust Bivariate Control Charts”**

### Authors

Vijayalakshmi S, Nicy Sebastian, Sajesh T.A

---

## Repository Contents

```text
Robust-Multivariate-Dispersion-Monitoring/
│
├── README.md
├── 01_MEWMS_UCL.R
├── 02_MEWMV_UCL.R
├── 03_MEWMS_Plot.R
├── 04_MEWMV_Plot.R
├── 05_MEWMV_ARL_Calculation.R
└── 06_MEWMS_ARL_Calculation.R
```

## R Scripts

### `01_MEWMS_UCL.R`

Contains the Monte Carlo simulation code used to obtain the empirical control limits for the MEWMS-based dispersion monitoring statistics.

### `02_MEWMV_UCL.R`

Contains the Monte Carlo simulation code used to obtain the control limits for the MEWMV-based dispersion monitoring statistics.

### `03_MEWMS_Plot.R`

Contains the code for calculating and plotting the MEWMS-based monitoring statistics for the covariance estimators considered in the study.

### `04_MEWMV_Plot.R`

Contains the code for calculating and plotting the MEWMV-based monitoring statistics for the covariance estimators considered in the study.

### `05_MEWMV_ARL_Calculation.R`

Contains the simulation code used to evaluate the Average Run Length (ARL) performance of the **MEWMV-based control charts** for the covariance estimators considered in the study.

The function

```r
BivIndCC_meanshift_combARL1()
```

implements the MEWMV-based ARL calculation.

### `06_MEWMS_ARL_Calculation.R`

Contains the simulation code used to evaluate the Average Run Length (ARL) performance of the **MEWMS-based control charts** for the covariance estimators considered in the study.

The function

```r
BivIndCCsim_combinedARL1()
```

implements the MEWMS-based ARL calculation.

---

## Data Source

The data used in this study are the **Chemical Process Data** reported in:

Montgomery, D. C. (2009). *Introduction to Statistical Quality Control*. 6th ed. John Wiley & Sons.

The Chemical Process Data are given in:

**Chapter 11: Multivariate Process Monitoring and Control, p. 520.**

The original data are not reproduced in this repository because they are available in the published reference cited above. Readers wishing to reproduce the analyses should obtain the Chemical Process Data from Montgomery (2009), Chapter 11, p. 520, and enter the data as specified in the relevant R scripts.

---

## Software Requirements

The analyses were performed using the **R statistical computing environment**.

The R scripts use standard R functions together with functions from packages including:

* `MASS`
* `robustbase`
* `rrcov`
* `pcaPP`
* `LaplacesDemon`
* `mvtnorm`
* `Matrix`
* `beepr`

The required packages can be installed in R using:

```r
install.packages(c(
  "MASS",
  "robustbase",
  "rrcov",
  "pcaPP",
  "LaplacesDemon",
  "mvtnorm",
  "Matrix",
  "beepr"
))
```

Package requirements may vary slightly depending on the particular script being executed.

---

## How to Reproduce the Results

1. Obtain the Chemical Process Data from Montgomery (2009), Chapter 11, p. 520.

2. Install R and the required R packages.

3. Download or clone this repository.

4. Enter the Chemical Process Data in the relevant R script as indicated in the code.

5. Run the UCL calculation scripts:

   * `01_MEWMS_UCL.R` for MEWMS control limits.
   * `02_MEWMV_UCL.R` for MEWMV control limits.

6. Use the resulting control limits in the corresponding plotting scripts:

   * `03_MEWMS_Plot.R`
   * `04_MEWMV_Plot.R`

7. Use the ARL simulation scripts to reproduce the ARL results reported in the manuscript:

   * `05_MEWMV_ARL_Calculation.R` for MEWMV-based ARL calculations.
   * `06_MEWMS_ARL_Calculation.R` for MEWMS-based ARL calculations.

8. Run the scripts sequentially as required for the particular analysis.

The scripts generate the numerical results and figures used in the manuscript.

---

## Reproducibility of Control Limits and ARL Results

The control limits and ARL values reported in the manuscript are obtained through simulation.

The control-limit scripts use Monte Carlo simulation to obtain the required empirical control limits. The ARL scripts use repeated simulation to evaluate the run-length performance of the proposed and competing control charts.

Because the results are simulation-based, small numerical differences may occur when the code is executed under different R versions, package versions, random-number-generation settings, or computational environments.

---

## Code Availability

The complete R code used for the analyses reported in the manuscript is provided in this repository to facilitate independent verification and reproduction of the results.

**GitHub repository:**
https://github.com/sajeshta/Robust-Multivariate-Dispersion-Monitoring

---

## Reproducibility Statement

The repository provides the complete computational code used for:

* estimation of control limits;
* computation of MEWMS monitoring statistics;
* computation of MEWMV monitoring statistics;
* generation of control-chart plots;
* simulation-based ARL evaluation; and
* reproduction of the numerical results reported in the manuscript.

The Chemical Process Data used in the study are identified by their published source so that readers can obtain the original data and reproduce the analyses using the accompanying R code.

---

## Reference

Montgomery, D. C. (2009). *Introduction to Statistical Quality Control*. 6th ed. John Wiley & Sons.
