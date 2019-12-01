## Local likelihood estimation of complex tail dependence structures, applied to U.S. precipitation extremes
## Daniela Castro-Camilo and Raphael Huser


We provide the R code needed to reproduce the simulation study. All results are obtained using R and the packages "mvtnorm", "fields", "numDeriv", "condMVNorm", "Matrix", "matrixcalc", "parallel", "plot3D", and "fda" all available on CRAN. Main code (with example usage) are denoted by [main]. We also provide a docker image to run the example usage in each main code.


### To fit the model

- [main] SimFit/SimFit.R: fit the model to the simulated data. The code runs on a single node in parallel across multiple cores. Includes an example usage and suggestions for extension to multiple nodes.
- [main] SimFit/TablesFigures.R: Examples to reproduce Tables and Figure in the simulation study.
- SimFit/CensoredLocalLikelihood.R: auxiliary code that contains the censored local log-likelihood functions.
- SimFit/logCopula.R, SimFit/logLikelihood.R, SimFit/logPartial.R: auxiliary code containing functions to compute the log copula for fully censored observations, the log-likelihood for non-censored observations, and the log partial derivatives of the copula for partially censored observations, respectively.
- SimFit/Tools.R: auxiliary functions.



### Docker image
We provide a docker image to run the example usage contained in each main code. For more info visit https://cloud.docker.com/u/daniccodes/repository/docker/daniccodes/local_lik_fcm_usprcpextremes. The code can be run, e.g., with:

- [To fit the model] `docker run --rm -ti danitest:latest Rscript SimulationStudy/SimFit/SimFit.R`
- [To print tables and produce figures] `docker run --rm -ti -v $PWD/Results/:/llFCM/Results/ danitest:latest Rscript SimulationStudy/SimFit/TablesFigures.R`




Daniela Castro-Camilo<br/>
School of Mathematics and Statistics<br/>
University of Glasgow<br/>
Glasgow G12 8QQ<br/>
UK

E-mail: daniela.castro.camilo@gmail.com




