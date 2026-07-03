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

source("functions.R")
source("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/Mod_biLocPol.R")


#### First run 'weather temperature/data_sets.R' ##### to obtain 'weather temperature/Berlin_lag0_setup.RData' ####
load("weather temperature/Berlin_lag0_setup.RData")                                                            ####
###################################################################################################################

p.eval = 96
x.eval.grid = 0:p.eval/p.eval


scale.x  = c(x.eval.grid[1] , x.eval.grid[41], x.eval.grid[81])
labels.x = c("00:00", "10:00","20:00")
scale.y  = c(x.eval.grid[41], x.eval.grid[81])
labels.y = c("10:00","20:00")


############################
###### marginal (lag 0)  ###
############################


x.eval.grid = tibble(TIME = hms::as_hms(c(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min"), as.POSIXct("1970-01-01 23:59:59"))))
quad.w      = simpson_weights(p.eval+1)

result.jan.est  = eigen( diag(sqrt(quad.w)) %*% jan.cov.lag0  %*% diag(sqrt(quad.w)), symmetric = TRUE) 
eig.val.jan.est = result.jan.est$values[1:3]
eig.vec.jan.est = result.jan.est$vectors[,1:3] * sqrt(quad.w) 

result.aug.est  = eigen( diag(sqrt(quad.w)) %*% aug.cov.lag0  %*% diag(sqrt(quad.w)), symmetric = TRUE) 
eig.val.aug.est = result.aug.est$values[1:3]
eig.vec.aug.est = result.aug.est$vectors[,1:3] * sqrt(quad.w) 


# calculate PCBF estimates for each month 
pcbf.jan.est    = data.frame(jan.cov.lag0   %*% ( eig.vec.jan.est[,1:3]  %*% diag(1/eig.val.jan.est) %*% diag(c(-1,-1, 1) ) )) 
pcbf.aug.est    = data.frame(aug.cov.lag0   %*% ( eig.vec.aug.est[,1:3]  %*% diag(1/eig.val.aug.est) %*% diag(c(-1, 1, 1) ) )) 



###############################################################################
###############################################################################
########                                                               ########
########                Dependent Multiplier Bootstrap                 ########
########                                                               ########
###############################################################################
###############################################################################
###############################################################################
start.24 = 21; end.24 = 117 # reduce to 24h data length (design grid 25) ######
N = 1000                    # Repetitions of procdure                    ######
###############################################################################

options(future.globals.maxSize = 64 * 1024^3) 
plan(multisession, workers = 30) 
  
# January
jan.mb = mb.PCBF(N, jan.data , list(w0), jan.cov.lag0, quadrature = "Simpson", k.max = 3, dependent = T, periodic = T, m = 144)
saveRDS(jan.mb, paste0("weather temperature/Results/jan.mb.rds") )
print("January done")
  
# August
aug.mb = mb.PCBF(N, aug.data , list(w0), aug.cov.lag0, quadrature = "Simpson", k.max = 3, dependent = T, periodic = T, m = 144)
saveRDS(aug.mb, paste0("weather temperature/Results/aug.mb.rds") )
print("August done")


plan(sequential)


##########################################################
#  Estimating quantiles and construct confindence bands  #
##########################################################

jan.crit.val = c()
aug.crit.val = c()

jan.mb.EW = readRDS(paste0("weather temperature/Results/jan.mb.rds"))$eig.val
aug.mb.EW = readRDS(paste0("weather temperature/Results/aug.mb.rds"))$eig.val

jan.sgn = c(-1,-1, 1) 
aug.sgn = c(-1, 1, 1) 

for(kk in 1:3){
  
  part = paste0("V",kk)
  
  n = dim(jan.data)[1]
  jan.mb.sample = readRDS(paste0("weather temperature/Results/jan.mb.rds"))$eig.vec[[part]]
  jan.mb.sd     = apply(jan.mb.sample,1,sd)
  jan.mb.q      = quantile( apply(abs(jan.mb.sample - pcbf.jan.est[,kk])/jan.mb.sd, 2, max), probs = 0.95, type = 2)
  jan.crit.val  = c(jan.crit.val,jan.mb.q * jan.mb.sd / sqrt(n))
  
  n = dim(aug.data)[1]
  aug.mb.sample = readRDS(paste0("weather temperature/Results/aug.mb.rds"))$eig.vec[[part]]
  aug.mb.sd     = apply(aug.mb.sample,1,sd)
  aug.mb.q      = quantile( apply(abs(aug.mb.sample - pcbf.aug.est[,kk])/aug.mb.sd, 2, max), probs = 0.95, type = 2)
  aug.crit.val  = c(aug.crit.val,aug.mb.q * aug.mb.sd / sqrt(n))

}


colnames(pcbf.jan.est)  = c("1", "2", "3")
colnames(pcbf.aug.est)  = c("1", "2", "3")

pcbf.jan.est    = pcbf.jan.est |> pivot_longer(cols = everything(), names_to = "PCBF",values_to = "value") |> arrange(PCBF)   |> mutate(eval.grid = rep(x.eval.grid$TIME, times = 3) )
pcbf.aug.est    = pcbf.aug.est |> pivot_longer(cols = everything(), names_to = "PCBF",values_to = "value") |> arrange(PCBF)   |> mutate(eval.grid = rep(x.eval.grid$TIME, times = 3) )


pcbf.jan.est = cbind(pcbf.jan.est, crit.val = jan.crit.val) |> mutate(UP = value + crit.val, LOW = value - crit.val)
pcbf.aug.est = cbind(pcbf.aug.est, crit.val = aug.crit.val) |> mutate(UP = value + crit.val, LOW = value - crit.val)



###############################################
# Plot of principal component basis functions #
###############################################

pcbf.est = rbind( data.frame(pcbf.jan.est, MONTH = factor("January", levels = c("January","August") )) , data.frame(pcbf.aug.est, MONTH = factor("August", levels = c("January","August") )) )

pcbf.colors = c("1" = "red2","2" = "blue3", "3" = "green4")


ggplot(data = pcbf.est) +
  
  geom_line(  aes(x = eval.grid, y = value, col = factor(PCBF), linetype = factor(PCBF)),   size = 1.2, show.legend = T, alpha = 1) +
  geom_ribbon(aes(x = eval.grid, ymin = LOW, ymax = UP, group = factor(PCBF)),  fill  = "grey6",  col = NA, alpha = 0.2,size = 0.5) +
  
  labs(x = " ", y = " ",title = bquote("p = 145 and  m = " *.(p.eval+1))) +
  scale_color_manual(values = pcbf.colors, name = "k") +
  scale_linetype_manual(values = c("dotdash", "dashed", "dotted"),name = "k") +
  
  theme(plot.title   = element_text(size = 22),
        legend.text  = element_text(size = 22),
        strip.text   = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.text.x  = element_text(size = 14),     
        axis.text.y  = element_text(size = 14),
        
        legend.key.size = unit(0.75, "cm")) +
  
  scale_y_continuous(breaks = c(-1.5, 0, 1.5), limits = c(-2, 2)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00")) +
  facet_grid(. ~ MONTH)

ggsave(paste0("Figures/estimation_marginal_PCBF.png"), width = 22, height = 14, units = "cm", dpi = 300)








####################
###### Long run  ###
####################

#### First run 'weather temperature/data_sets.R' ##### to obtain 'weather temperature/Berlin_lagb_setup.RData' ####
load("weather temperature/Berlin_lagb_setup.RData")                                                            ####
###################################################################################################################


x.eval.grid = tibble(TIME = hms::as_hms(c(seq(from = as.POSIXct("1970-01-01 00:00:00"),to   = as.POSIXct("1970-01-01 23:45:00"),by   = "15 min"), as.POSIXct("1970-01-01 23:59:59"))))
quad.w      = simpson_weights(p.eval+1)

result.jan.est  = eigen( diag(sqrt(quad.w)) %*% lr.Gamma.d.month[[1]]  %*% diag(sqrt(quad.w)), symmetric = TRUE) 
eig.val.jan.est = result.jan.est$values[1:3]
eig.vec.jan.est = result.jan.est$vectors[,1:3] * sqrt(quad.w) 

result.aug.est  = eigen( diag(sqrt(quad.w)) %*% lr.Gamma.d.month[[8]]  %*% diag(sqrt(quad.w)), symmetric = TRUE) 
eig.val.aug.est = result.aug.est$values[1:3]
eig.vec.aug.est = result.aug.est$vectors[,1:3] * sqrt(quad.w) 


# calculate PCBF estimates for each month 
pcbf.jan.est    = data.frame(lr.Gamma.d.month[[1]]   %*% ( eig.vec.jan.est[,1:3]  %*% diag(1/eig.val.jan.est) %*% diag(c(-1,1, 1))) ) 
pcbf.aug.est    = data.frame(lr.Gamma.d.month[[8]]   %*% ( eig.vec.aug.est[,1:3]  %*% diag(1/eig.val.aug.est) %*% diag(c(-1,1,-1))) ) 




#############################################################################
# Staring multiplier botstrap with dependent multipliers ####################
#############################################################################
N = 1000                    # Repetitions of procdure                    ####
#############################################################################

options(future.globals.maxSize = 64 * 1024^3) 
plan(multisession, workers = 30) 

# January
jan.mb.LR = mb.PCBF(N, jan.data , list(w0, w1, w2), lr.Gamma.d.month[[1]], quadrature = "Simpson", k.max = 3, dependent = T, periodic = T, m = 144, max.lag = max.lag[1,2])
saveRDS(jan.mb.LR, paste0("weather temperature/Results/jan.mb.LR.rds") )
print("January done: long run kernel")

# August
aug.mb.LR = mb.PCBF(N, aug.data , list(w0, w1, w2), lr.Gamma.d.month[[8]], quadrature = "Simpson", k.max = 3, dependent = T, periodic = T, m = 144, max.lag = max.lag[8,2])
saveRDS(aug.mb.LR, paste0("weather temperature/Results/aug.mb.LR.rds") )
print("August done: long run kernel")

plan(sequential)


##########################################################
#  Estimating quantiles and construct confindence bands  #
##########################################################

jan.crit.val = c()
aug.crit.val = c()

jan.mb.EW = readRDS(paste0("weather temperature/Results/jan.mb.LR.rds"))$eig.val
aug.mb.EW = readRDS(paste0("weather temperature/Results/aug.mb.LR.rds"))$eig.val

jan.sgn = c(-1, 1, 1) 
aug.sgn = c(-1, 1,-1) 

for(kk in 1:3){
  
  part = paste0("V",kk)
  
  n = dim(jan.data)[1]
  jan.mb.sample = readRDS(paste0("weather temperature/Results/jan.mb.LR.rds"))$eig.vec[[part]]
  jan.mb.sd     = apply(jan.mb.sample,1,sd)
  jan.mb.q      = quantile( apply(abs(jan.mb.sample - pcbf.jan.est[,kk])/jan.mb.sd, 2, max), probs = 0.95, type = 2)
  jan.crit.val  = c(jan.crit.val,jan.mb.q * jan.mb.sd / sqrt(n))
  
  n = dim(aug.data)[1]
  aug.mb.sample = readRDS(paste0("weather temperature/Results/aug.mb.LR.rds"))$eig.vec[[part]]
  aug.mb.sd     = apply(aug.mb.sample,1,sd)
  aug.mb.q      = quantile( apply(abs(aug.mb.sample - pcbf.aug.est[,kk])/aug.mb.sd, 2, max), probs = 0.95, type = 2)
  aug.crit.val  = c(aug.crit.val,aug.mb.q * aug.mb.sd / sqrt(n))
  
}


colnames(pcbf.jan.est)  = c("1", "2", "3")
colnames(pcbf.aug.est)  = c("1", "2", "3")

pcbf.jan.est    = pcbf.jan.est |> pivot_longer(cols = everything(), names_to = "PCBF",values_to = "value") |> arrange(PCBF)   |> mutate(eval.grid = rep(x.eval.grid$TIME, times = 3) )
pcbf.aug.est    = pcbf.aug.est |> pivot_longer(cols = everything(), names_to = "PCBF",values_to = "value") |> arrange(PCBF)   |> mutate(eval.grid = rep(x.eval.grid$TIME, times = 3) )

pcbf.jan.est = cbind(pcbf.jan.est, crit.val = jan.crit.val) |> mutate(UP = value + crit.val, LOW = value - crit.val)
pcbf.aug.est = cbind(pcbf.aug.est, crit.val = aug.crit.val) |> mutate(UP = value + crit.val, LOW = value - crit.val)



###############################################
# Plot of principal component basis functions #
###############################################

pcbf.est = rbind( data.frame(pcbf.jan.est, MONTH = factor("January", levels = c("January","August") )) , data.frame(pcbf.aug.est, MONTH = factor("August", levels = c("January","August") )) )


ggplot() +
  
  geom_line(mapping = aes(x = eval.grid, y = value, col = factor(PCBF), linetype = factor(PCBF)), data = pcbf.est,   size = 1.2, show.legend = T, alpha = 1) +
  geom_ribbon(aes(x = eval.grid, ymin = LOW, ymax = UP, group = factor(PCBF)), data = pcbf.est,  fill  = "grey6",  col = NA, alpha = 0.2,size = 0.5) +
  
  labs(x = " ", y = " ",title = bquote("p = 145 and  m = " *.(p.eval+1))) +
  scale_color_manual(values = pcbf.colors, name = "k") +
  scale_linetype_manual(values = c("dotdash", "dashed", "dotted"),name = "k") +
  
  theme(plot.title   = element_text(size = 22),
        legend.text  = element_text(size = 22),
        strip.text   = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.text.x  = element_text(size = 14),     
        axis.text.y  = element_text(size = 14),
        
        legend.key.size = unit(0.75, "cm")) +
  
  scale_y_continuous(breaks = c(-1.5, 0, 1.5), limits = c(-2, 2)) +
  scale_x_time(breaks = as_hms(c("00:00:00", "10:00:00", "20:00:00")),labels = c("00:00", "10:00", "20:00")) +
  facet_grid(. ~ MONTH)

ggsave(paste0("Figures/estimation_long_run_PCBF.png"), width = 22, height = 14, units = "cm", dpi = 300)









