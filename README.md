# BSOINN
 * [BSOINN-package](#BSOINN-package)
    1. [Overview](#Overview)
    2. [Installing BSOINN](#InstallingBSOINN)
    3. [Democode-simulation](#Democode-simulation)
    4. [Democode-realdata](#Democode-realdata)
    5. [Implementing Democode](#implementingDemocode)

# BSOINN-package <a name="BSOINN-package"></a>

"BSOINN" package conducts Bayesian infernece for a scalar on image regression model with non-ignorable non-response, which implements the method proposed in the paper "Bayesian scalar on image regression with non-ignorable non-response".

## Overview  <a name="Overview"></a>

This package contains three functions for the Bayesian sclar on image (BSOI) regression model. 

The three funcitons of the package are as follows:

1. The function "BSOIFull" conducts Bayesian inference for a BSOI with fully-observed response.

2. The function "BSOIIN" conducts Bayesian inference for a BSOI with ignorable non-response.

3. The function "BSOINN" conducts Bayesian inference for a BSOI with non-ignorable non-response.

The three functions are based on the funcitonal principal component analysis (FPCA) for imaging data, which can be implemented using the "fast.svd" function in the package of "corpcor". Details of the FPCA can be found in the paper and the democode for the simulation study. 

##Intalling BSOINN  <a name="InstallingBSOINN"></a>

* Installing BSOINN on Windows:
1. Download and install R software from http://cran.r-project.org.
2. Download and install Rtools from  http://cran.r-project.org/bin/windows/Rtools/. During the installation process, please check the box provided to edit the system PATH, so that the C++ compliler included in Rtools can be used by R.
3. Download and install Rstudio software from https://www.rstudio.com/.
4. Install packages "Rcpp" and "RcppArmadillo" from CRAN inside the R software.
5. Install package "BSOINN" from the local package archive "BSOINN_1.0.tar.gz".

* Installing BSOINN on Mac OS:
1. Download and install R software.
2. Download and install C++ Toolchain, Xcode, from Mac "App Store". After the Xcode is installed, you need to open Xcode once to accept the license aggrement.
3. Download and install Rstudio software.
4. Run the following code in R to put the path of "libgfortran"" into "FLIBS": 
```{r, eval=FALSE, warning=FALSE,message=FALSE}
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("FLIBS = -L`gfortran -print-file-name=libgfortran.dylib | xargs dirname`", file = M, sep = "\n", append = TRUE)
```
5. Install packages "Rcpp" and "RcppArmadillo" from CRAN inside the R software.
6. Install package "BSOINN" from the local package archive "BSOINN_1.0.tar.gz".

##Democode-simulation <a name="Democode-simulation"></a>

The democode for the simulation study from the paper can be found in the file "simulation_example.R" from the folder "Democode-simulation". 

After installing the package "BSOINN" (this package) and "corpcor" (can be installed from CRAN), this R program can be implemented. 

Please refer to the R file for detailed codes and explanations.

## Democode-realdata <a name="Democode-realdata"></a>

The democode for the real data analysis using ADNI data set from the paper can be found in the file "realdata_example.R" from the folder "Democode-realdata". 

Data from the ADNI study are inluded:

* "Y.txt": learning scores of the subjects at the 36th month, with non-response.

* "COV.txt": gender, age, education (the instrument),  race, whether ever married, apoe gene, whether MCI or AD.

* "XI.txt": five eigenscores corresponding to the first five eigenimages of the imaging data.

* "R.txt": 0 denotes response and 1 denotes non-response.

After installing the package "BSOINN" (this package), this R program can be implemented. 

Please refer to the R file for detailed codes and explanations.

## Implementing Democode <a name="implementingDemocode"></a>

* Democode-simulation:

1. set the folder "Democode-simulation" as the working directory in R.

2. implementing the file "simulation_example.R" through the following command:
```{r, eval=FALSE, warning=FALSE,message=FALSE}
source("simulation_example.R")
```

3. The convergence of the algorithm is demonstrated through the trace-plot of three parallel chains, and the estimation results are automatically shown on the R console.

* Democode-realdata:

1. set the folder "Democode-realdata" as the working directory in R.

2. implementing the file "realdata_example.R" through the following command:
```{r, eval=FALSE, warning=FALSE,message=FALSE}
source("realdata_example.R")
```

3. The convergence of the algorithm is demonstrated through the trace-plot of three parallel chains, and the estimation results are automatically shown on the R console.


