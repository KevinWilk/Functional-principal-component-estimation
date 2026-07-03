##################################################################################
#############  Simulation  #######################################################
##################################################################################

library(MASS)
library(ggplot2)
library(locpol)
library(interp)
library(stats)
library(future.apply)
library(parallel)
library(biLocPol)
library(plotly)
library(tidyr)
library(crayon)
library(future)
library(hms)
library(dplyr)
library(tibble)
library(gridExtra)

source("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/Mod_biLocPol.R")
source("functions.R")

###############################################################
#k = 2 # For k = 1,2: functions of analytical solution ######
###############################################################

n.seq  = c(50,100,200,350,500)
p.seq  = c(15,25,50,75,100,125,150)
p.eval = 200 #(for Simpson rule)
n      = 500
N      = 1000 # number of repitions  

# Parameter OU Process
theta = 3; sigma = 2

# Standard deviation for additional errors
sd = 1

x.eval.grid = 0:p.eval/p.eval

####### bandwidth selection######################################################
H = lapply(1:length(p.seq), function(l){seq(1, max(0.05,3/p.seq[l]), -0.05)})####
#################################################################################






options(future.globals.maxSize = 64 * 1024^3) 
plan(multisession, workers = 60)               


#####################
# eigenvector optim #
#####################
for(i in 1:5){
for(j in 1:5){ 
    
    ####################
    # Trapezoidal rule #
    ####################
    
    # phi_1
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 1, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "trapezoidal")  

    # phi_2
    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 2, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "trapezoidal")   
    
    optim = list(eigfunc_1 = optim1, eigfunc_2 = optim2)
    saveRDS(optim, paste0("Results/eigenvector/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )
    
    print(paste0("Trapezoidal rule done: p = ",p.seq[i], " and n = ",n.seq[j]))
    

    ################
    # Simpson rule #
    ################

    # phi_1
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 1, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "Simpson")  

    # phi_2
    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 2, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "Simpson")    
    
    optim = list(eigfunc_1 = optim1, eigfunc_2 = optim2)
    saveRDS(optim, paste0("Results/eigenvector/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )

    print(paste0("Simpson rule done: p = ",p.seq[i], " and n = ",n.seq[j]))
    
}
}


####################
# covariance optim #
####################
for(i in 1:5){ 
    

    optim = h.optim.cov(N, n.seq[5], p.seq[i], p.eval, H[[i]], 
                        m = 1, w.parallel = T, 
                        f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd))  

    print(paste0("done: p = ",p.seq[i], " and n = ",n.seq[5]))

    saveRDS(optim, paste0("Results/covariance/Cov_p_",p.seq[i],"_n_",n.seq[5],"_optim.rds") )

}

############################
# supnorm error eigenvalue #
############################
for(i in 1:5){
for(j in 1:5){    
    
    ####################
    # Trapezoidal rule #
    ####################
    
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 1, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "trapezoidal", supnorm = "eigenvalue")  

    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 2, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "trapezoidal", supnorm = "eigenvalue")  

    optim = list(eigval_1 = optim1, eigval_2 = optim2)
    
    saveRDS(optim, paste0("Results/eigenvalue/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )
    print(paste0("Trapezoidal rule done: p = ", p.seq[i], " and n = ", n.seq[j]))


    
    ################
    # Simpson rule #
    ################
    
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 1, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "Simpson", supnorm = "eigenvalue")  

    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 2, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "Simpson", supnorm = "eigenvalue")  

    optim = list(eigval_1 = optim1, eigval_2 = optim2)
    
    saveRDS(optim, paste0("Results/eigenvalue/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )
    print(paste0("Simpson rule done: p = ", p.seq[i], " and n = ", n.seq[j]))
}
}








options(future.globals.maxSize = 64 * 1024^3) 
plan(multisession, workers = 20)               # p > 100 only 20 workers


#####################
# eigenvector optim #
#####################
for(i in 6:7){
for(j in 1:5){ 
    
    ####################
    # Trapezoidal rule #
    ####################
    
    # phi_1
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 1, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "trapezoidal")  

    # phi_2
    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 2, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "trapezoidal")   
    
    optim = list(eigfunc_1 = optim1, eigfunc_2 = optim2)
    saveRDS(optim, paste0("Results/eigenvector/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )
    
    print(paste0("Trapezoidal rule done: p = ",p.seq[i], " and n = ",n.seq[j]))
    

    ################
    # Simpson rule #
    ################

    # phi_1
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 1, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "Simpson")  

    # phi_2
    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                    k = 2, m = 1, w.parallel = T, 
                    f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                    quadrature = "Simpson")    
    
    optim = list(eigfunc_1 = optim1, eigfunc_2 = optim2)
    saveRDS(optim, paste0("Results/eigenvector/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )

    print(paste0("Simpson rule done: p = ",p.seq[i], " and n = ",n.seq[j]))
    
}
}



####################
# covariance optim #
####################
for(i in 6:7){ 
    

    optim = h.optim.cov(N, n.seq[5], p.seq[i], p.eval, H[[i]], 
                        m = 1, w.parallel = T, 
                        f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd))  

    print(paste0("done: p = ",p.seq[i], " and n = ",n.seq[5]))

    saveRDS(optim, paste0("Results/covariance/Cov_p_",p.seq[i],"_n_",n.seq[5],"_optim.rds") )

}


############################
# supnorm error eigenvalue #
############################
for(i in 6:7){
for(j in 1:5){    
    
    ####################
    # Trapezoidal rule #
    ####################
    
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 1, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "trapezoidal", supnorm = "eigenvalue")  

    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 2, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "trapezoidal", supnorm = "eigenvalue")  

    optim = list(eigval_1 = optim1, eigval_2 = optim2)
    
    saveRDS(optim, paste0("Results/eigenvalue/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )
    print(paste0("Trapezoidal rule done: p = ", p.seq[i], " and n = ", n.seq[j]))


    
    ################
    # Simpson rule #
    ################
    
    optim1 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 1, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "Simpson", supnorm = "eigenvalue")  

    optim2 = h.optim(N, n.seq[j], p.seq[i], p.eval, H[[i]], 
                     k = 2, m = 1, w.parallel = T, 
                     f = mu0, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd),
                     quadrature = "Simpson", supnorm = "eigenvalue")  

    optim = list(eigval_1 = optim1, eigval_2 = optim2)
    
    saveRDS(optim, paste0("Results/eigenvalue/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds") )
    print(paste0("Simpson rule done: p = ", p.seq[i], " and n = ", n.seq[j]))
}
}

plan(sequential)






