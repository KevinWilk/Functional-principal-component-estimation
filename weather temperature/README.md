## Overview 

#### `weather temperature.R` 
- runs dependent multiplier bootstrap on mariginal lag 0 and long-run kernel with `B=1000 repetitions` <br>
-> necessary data for mariginal lag 0: `weather temperature/Berlin_lag0_setup.RData`) and saves *L_2 normed* multiplier bootstrap samples and eigenvalues as `Results/month.mb.rds`  <br>
-> necessary data for long-run kernel (`weather temperature/Berlin_lagb_setup.RData`) and saves *L_2 normed* multiplier bootstrap samples and eigenvalues as `Results/month.mb.LR.rds` <br>
- constructs 95% uniform confidence bands <br>
- creates Figures 4 and 5 which are saved in *Figures/* <br>


#### `data_sets.R`
- creates necessary data for dependent multiplier bootstrap on mariginal lag 0 (`weather temperature/Berlin_lag0_setup.RData`) and long-run kernel (`weather temperature/Berlin_lagb_setup.RData`)


  

