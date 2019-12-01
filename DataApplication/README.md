## Local likelihood estimation of complex tail dependence structures, applied to U.S. precipitation extremes
## Daniela Castro-Camilo and Raphael Huser


We provide the R code needed to reproduce the application to U.S. precipitation extremes. All results are obtained using R and the packages "mvtnorm", "R.utils", "fields", "numDeriv", "condMVNorm", "Matrix", "matrixcalc", "rje", "ismev", "evd", and "parallel", "homtest", "ggmap", "ggplot2", and "reshape2" all available on CRAN. Main code (with example usage) are denoted by [main]. We also provide a docker image to run the example usage in each main code.


### To fit the model

- [main] Fit/FitUSprcp.R: fit our model to the U.S. precipitation extreme observations. The code runs on a single node in parallel across multiple cores. Includes an example usage and suggestions for extension to multiple nodes.
- Fit/CensoredLocalLikelihood.R: auxiliary code that contains the censored local log-likelihood functions.
- Fit/logCopula.R, Fit/logLikelihood.R, Fit/logPartial.R: auxiliary code containing functions to compute the log copula for fully censored observations, the log-likelihood for non-censored observations, and the log partial derivatives of the copula for partially censored observations, respectively.
- Fit/Tools.R: auxiliary functions.

- Data/grid60.txt: regular grid with 2235 grid points at an internodes distance of 60km, where the model will be fitted.
- Data/locations1218.txt: all the 1218 station's locations.
- Data/StartingValues.txt: starting values for the optimizer
- Data/datamat5days1218.Rdata: five days-cumulative precipitation in winter time. This is a 2070x1218 matrix.
- Data/neigh_grid60km_dist150km_forward.Rdata: neighbors id for each grid point obtained by testing homogeneity of the marginal distributions.



### To estimate uncertainty via Bootstrap

- Fit/Bootstrap/BlockBootstrap.R: generates block-bootstrap samples.
- [main] Fit/Bootstrap/BFitUSprcp.R: fit our model to the U.S. precipitation extreme observations, for every block-Bootstrap sample b.



### To generate figures 

Before running, introduce your api key in the <key> line to enable Google services in R. If using the docker image, then do Rscript Plot/EmpiricalQ.R "api_key<-'your_api_key'". Same for Plots/Chi.R and Plots/BoxplotsChisd.R.

- [main] Plots/EmpiricalQ.R: plot empirical quantiles of five day-cumulative winter precipitation data (Figure 6 in the paper).
- [main] Plots/Chi.R: plot estimates of $\chi_h(u)$ (Figure 7 in the paper).
- [main] Plots/BoxplotsChisd.R: boxplots of standard deviations for estimates of $\chi_h(u)$ (Figure 8 in the paper). Same code can be used to generate boxplots of standard deviations for estimates of log-rate and log-range (Figure 9 in the paper).
- Plots/param_nu=0.5.txt: estimates of the rate and range parameters obtained using Fit/FitUSprcp.R. Useful to run Plots/Chi.R.
- Plots/sdchi_nu=0.5.txt: estimates of sd of $\chi_h(u)$ for h = 20, 40, 80, 160 and u = 0.95. Useful to run Plots/BoxplotsChisd.R.



### To compute return periods (Table 3 in the paper)

- [main] ReturnPeriods/Parametric/ParametricReturnPeriods.R: Compute parametric return periods associated with catastrophic events over selected stations in 5 states. The return periods are computed over simulations from our model, and the code includes functions to simulate N replicates from our model and transformation to uniform margins.
- [main] ReturnPeriods/Parametric/BootstrapParametricReturnPeriods.R: same as ReturnPeriods/Parametric/ParametricReturnPeriods.R, for every Bootstrap sample. Includes Bootstrap-based standard deviations for parametric return period estimates.
- [main] ReturnPeriods/Nonparametric/NonParamReturnPeriods.R: Compute non-parametric return periods associated with catastrophic events over selected stations in 5 states. Includes Bootstrap-Based standard deviations.



### To obtain neighbors for each grid point (no need to run, neighbors are already computed. Included here for reproducibility purposes)

- [main] Neighbors/NeighborhoodSelection.R: main code to generated a homogeneous neighborhood for each grid point.
- Neighbors/myHtests.R Neighbors/HoTests.R: auxiliary functions.



### To preprocess the data (no need to run, data are already preprocessed. Included here for reproducibility purposes)

- [main] PreProcessing/1_PreProcessing.R: Code to remove observations with quality flags (run this first)
- [main] PreProcessing/2_5DaysCumRainfall.R: Code to put 5-days cumulative winter precipitation data in a N x n matrix (run this after PreProcessing/1_PreProcessing.R)
- Preprocessing/RawData: folder with the raw data (zip files) for each state. Need to unzip to run PreProcessing/1_PreProcessing.R.
- Preprocessing/States: Folder where the processed data obtained when running PreProcessing/1_PreProcessing.R will be saved.


###  Docker image
We provide a docker image to run the example usage contained in each main code. For more info visit https://cloud.docker.com/u/daniccodes/repository/docker/daniccodes/local_lik_fcm_usprcpextremes. The code can be run, e.g., with:

- [To fit the model] `docker run --rm -ti daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Fit/FitUSprcp.R`
- [To estimate uncertainty via Bootstrap] `docker run --rm -ti daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Fit/Bootstrap/BFitUSprcp.R`
- [To generate block-bootstrap samples] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Fit/Bootstrap/BlockBootstrap.R`
- [To generate figures] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Plots/EmpiricalQ.R "api_key<-'your_api_key'"`
- [To generate figures] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Plots/Chi.R "api_key<-'your_api_key'"`
- [To generate figures] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Plots/BoxplotsChisd.R "api_key<-'your_api_key'"`
- [To compute return periods] `docker run --rm -ti daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/ReturnPeriods/Nonparametric/NonParamReturnPeriods.R`
- [To obtain neighbors for each grid point. Takes time!] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/Neighbors/NeighborhoodSelection.R`
- [To preprocess the data. Please read the warning message!] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/PreProcessing/1_PreProcessing.R`
- [To preprocess the data. Please read the warning message!] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ daniccodes/local_lik_fcm_usprcpextremes:latest Rscript DataApplication/PreProcessing/2_5DaysCumRainfall.R`




Daniela Castro-Camilo<br/>
School of Mathematics and Statistics<br/>
University of Glasgow<br/>
Glasgow G12 8QQ<br/>
UK

E-mail: daniela.castro.camilo@gmail.com




