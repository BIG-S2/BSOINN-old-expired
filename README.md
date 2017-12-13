# BSOINN
Bayesian Scalar on Image Regression with Non-ignorable Non-response

Main functions are included in the R package "BSOINN", please install the package first to implement the demo codes.

Please note that the brain images are too large to be included in the package or in the folder of demo code for the real data analysis. We therefore require the eigenscores from the functional principal component analysis (FPCA) on the images as inputs for the functions in the “BSOINN” package. The FPCA in the paper can be easily implemented using matlab or the “fast.svd” function in the R, which is demonstrated in the demo code for the simulation study. 


Please also install the package "corpcor" from CRAN to use the function "fast.svd"

Demo code for simulation is in the folder "Democode-simulation".

The simulation_example.R generate simulated data sets and conduct the analysis.



Demo code the real data analysis is in the folder "Democode-realdata".

The realdata_example.R read data (based on the ADNI study) from files :

--- Y.txt: learning scores of the subjects at the 36th month, with missingness

--- COV.txt:  gender, age, education, race, whether ever married, apoe4,  whether MCI or AD  

--- XI.txt: five eigenscores corresponding to the first five eigenimages of RAVENS data.

--- R.txt: 0, observed, 1, missing
