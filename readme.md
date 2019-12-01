# Local Likelihood Estimation of Complex Tail Dependence Structures, Applied to U.S. Precipitation Extremes

# Author Contributions Checklist Form

## Data

### Abstract

Database defined in a network of 1218 stations in the contiguous US comprising key variables
for climate monitoring, such as precipitation, snowfall, maximum and minimum temperature,
etc. The observations were daily recorded from the beginning of the twentieth century till 2014.

### Availability 

The raw data are publicly available from the USHCN website, but we also made them public in a
GitHub repository, along with the cleaned data. Please visit https://github.com/DaniCCodes/
LocalLikFCM_USprcpExtremes

### Description 

Permissions (demonstrate that author has legitimate access to data): the data are public.

Licensing information: not applicable, the data are public.

Link to data: http://cdiac.ess-dive.lbl.gov/ftp/ushcn_daily/

Data provenance, including identifier or link to original data if different than above: US
historical climatology archive (USHCN, same link as above).

File format: .txt.gz files.

Metadata: http://cdiac.ess-dive.lbl.gov/epubs/ndp/ushcn/daily_doc.html

Version information: not available.


## Code

### Abstract (Mandatory)

The repository contains the required R codes to replicate the simulation study, the data
application, figures, and tables in "Local likelihood estimation of complex tail dependence
structures, applied to U.S. precipitation extremes."

### Description

How delivered (R package, Shiny app, etc.): R codes in a public repository.

Licensing information (default is MIT License): MIT license.

Link to code/repository: https://github.com/DaniCCodes/LocalLikFCM_USprcpExtremes

Version information (e.g., for a Git repository, the number or branch+commit): GitHub
repository 1 branch, 23 commits.

### Optional Information

Hardware requirements (e.g., operating system with version number, access to cluster, GPUs,
etc.): access to cluster is advisable.

Supporting software requirements (e.g., libraries and dependencies, including version numbers):

R version 3.5.2. 

Libraries and their versions in parenthesis: mvtnorm (1.0-10), R.utils (2.8.0),
fields (9.8-1), numDeriv (2016.8-1), condMVNorm (2015.2-1), Matrix (1.2-17), matrixcalc
(1.0-3), rje (1.9), ismev (1.42), evd (2.3-3), parallel (base), homtest (1.0-5), ggmap (3.0.0),
ggplot2 (3.1.1), and reshape2 (1.4.3).

## Instructions for Use

### Reproducibility 

What is to be reproduced: the simulation study, the data application, figures, and tables in
"Local likelihood estimation of complex tail dependence structures, applied to U.S.
precipitation extremes."

How to reproduce analyses: the set of codes include example usage.
