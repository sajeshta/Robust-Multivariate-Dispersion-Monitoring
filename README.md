**Reproducibility Materials**

This repository contains the complete R code used to reproduce the analyses, numerical results, tables, and figures reported in the manuscript.

**Manuscript**

Title: Monitoring Process Variability with Robust Bivariate Control Charts

Authors: Vijayalakshmi S, Nicy Sebastian, Sajesh T.A

**Repository Contents**

├── README.md

└── R/

    └── analysis.R
    
R/

This folder contains the R code used for data preparation, statistical analysis, simulation studies, computation of performance measures, and generation of the results reported in the manuscript.

**Data Source**

The data used in this study are the Chemical Process Data reported in:

Montgomery, D. C. (2009). Introduction to Statistical Quality Control. 6th ed. John Wiley & Sons.

The Chemical Process Data are given in Chapter 11: Multivariate Process Monitoring and Control, p. 520.

The original data are not reproduced in this repository because they are available in the published reference cited above.
Readers wishing to reproduce the analyses should obtain the Chemical Process Data from Montgomery (2009), Chapter 11, p. 520, and use the data as specified in the R code provided in this repository.

**Software Requirements**

The analysis was performed using the R statistical computing environment.

The R script specifies the packages required for the analysis. Required packages can be installed using:

install.packages("package_name")

**How to Reproduce the Results**

1.	Obtain the Chemical Process Data from Montgomery (2009), Chapter 11, p. 520.
2.	Install R and the required R packages.
3.	Download or clone this repository.
4.	Open the R script located in the R/ folder.
5.	Enter the Chemical Process Data as indicated in the R script.
6.	Run the R script sequentially.
7.	The code produces the numerical results, tables, and figures reported in the manuscript.
   
**Code Availability**

The complete R code used to perform the analyses reported in the manuscript is provided in this repository to facilitate independent verification and reproduction of the results.

**Reproducibility**

The code has been organized to allow the analyses reported in the manuscript to be reproduced using the Chemical Process Data obtained from Montgomery (2009).

Minor differences in numerical results may occur because of differences in R versions, package versions, or computational environments.

**Reference**

Montgomery, D. C. (2009). Introduction to Statistical Quality Control. 6th ed. John Wiley & Sons.

