## Overview: Simulation 

#### `supnorm_error.R` with 1000 repetions: calculation of sup-norm error of PCBFs (k = 1,2) and covariance kernel and in addition of error of corresponding eigenvalues
- `n.seq  = c(50,100,200,350,500)` and `p.seq  = c(15,25,50,75,100,125,150)` <br>
- `optim`: as list with elements optim1 (k = 1: list of sup-norm error for each h) and optim2 (k = 2: list of sup-norm error for each h) <br>
  -> saved as `Trapez/Simpson_p_p.seq[i]_n_n.seq[j]_optim.rds"` (trapezodial or Simpson rule) in *Results/...* <br>
  -> *.../eigenvector/* for PCBFs, *.../eigenvalue/* for corresponding eigenvalues or *.../covariance/* for covariance kernel


  

