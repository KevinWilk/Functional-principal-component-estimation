##################################################################################
#############  Figures  ##########################################################
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







#########################################################################
#########################################################################
#####                                  ##################################
#####  Illustrativ estimation of PCBF  ##################################
#####                                  ##################################
#########################################################################
#########################################################################

load("Figures/illustrativ_example.RData")

# Parameters of design.grid and observational rows
p.seq  = c(15,25,50,75,100,145)
n      = 500
i      = 3 # p = 50

# Quadrature points
p.eval      = 200 #(for Simpson rule)
x.eval.grid = 0:p.eval/p.eval


k = 2 # For k = 1,2: functions of analytical solution 

# Parameter OU Process
theta = 3; sigma = 2

# determine the roots of function `cond()`
omega_k_list  = sapply(1:k, omega_k, theta = theta)

# determine the corresponding analytical eigenvalues lambda_k 
lambda_k_list = sapply(omega_k_list, function(x, sigma, theta){sigma^2/(theta^2+x^2)}, sigma = sigma, theta = theta)

# determine the corresponding analytical principal component basis function phi_k
phi_k_list    = sapply(omega_k_list, function(x, omega, theta){ phi_k(x , omega, theta)  }, theta = theta, x = x.eval.grid) |>
  as.data.frame() |>
  #mutate(across(everything(), ~ .x / sqrt(sum(.x^2)))) |> #scale to discretisation
  mutate(x = x.eval.grid) |>
  setNames(c("1", "2","eval.grid")) |>
  pivot_longer(cols = c(starts_with("1"), starts_with("2")),names_to = "PCBF", values_to = "value") |>
  arrange(PCBF)



# Standard deviation for additional errors
sd = 1

# Gernate data set
x.design = 0:p.seq[i]/p.seq[i]
set.seed(2312)
Y = FDA_observation(n, x.design, f = mu, r.process = OU, process.arg = list(alpha = theta, sigma = sigma), eps.arg = list(mean = 0, sd = sd))

###############################################
# estimate pricipal compenent basis functions #
###############################################

# compute bivariate kernel weights
options(future.globals.maxSize = 20 * 1024^3) 
plan(multisession, workers = 20) 

bw.1 = 0.65 # optimal for Simpson rule
w.1  = local.polynomial.weights(p.seq[i]+1, bw.1, p.eval = p.eval+1, parallel = T, m = 1, del = 0, grid.type = "less",   eval.type = "full", parallel.environment = F) # from Wilk and Holzmann (2026)

bw.2 = 0.25 # optimal for Simpson rule
w.2  = local.polynomial.weights(p.seq[i]+1, bw.2, p.eval = p.eval+1, parallel = T, m = 1, del = 0, grid.type = "less",   eval.type = "full", parallel.environment = F) # from Wilk and Holzmann (2026)

plan(sequential)

# calculate convariance matrix estimation 
cov.est.1 = eval.weights(w.1, observation.transformation(Y, grid.type = "less")) # from Wilk and Holzmann (2026)
cov.est.2 = eval.weights(w.2, observation.transformation(Y, grid.type = "less")) # from Wilk and Holzmann (2026)

#quad.w  = trapezoid_weights(p.eval+1)
quad.w  = simpson_weights(p.eval+1)

result1.est  = eigen( diag(sqrt(quad.w)) %*% cov.est.1  %*% diag(sqrt(quad.w)), symmetric = TRUE) 
eig.val1.est = result1.est$values[1]
eig.vec1.est = result1.est$vectors[,1] * sqrt(quad.w) # quad.w (sqrt(quad.w) * v )  ,see Kress (2014) equation (12.9)
 
result2.est  = eigen( diag(sqrt(quad.w)) %*% cov.est.2  %*% diag(sqrt(quad.w)), symmetric = TRUE) 
eig.val2.est = result2.est$values[2]
eig.vec2.est = result2.est$vectors[,2] * sqrt(quad.w) # quad.w (sqrt(quad.w) * v )  ,see Kress (2014) equation (12.9)



# calculate for each estimator of PCBF its sign with analytical PCBF
sign.est    = diag(sign(crossprod(sapply(omega_k_list, function(x, omega, theta){ phi_k(x , omega, theta)  }, theta = theta, x = x.eval.grid), cbind(eig.vec1.est,eig.vec2.est))))
pcbf.est    = data.frame(cov.est.1   %*% ( eig.vec1.est  * 1/eig.val1.est * sign.est[1]), cov.est.2   %*% ( eig.vec2.est  * 1/eig.val2.est * sign.est[2])) 
colnames(pcbf.est)  = c("1", "2")
pcbf.est    = pcbf.est |> pivot_longer(cols = everything(), names_to = "PCBF",values_to = "value") |> arrange(PCBF)   |> mutate(eval.grid = rep(x.eval.grid, times = 2) )
pcbf.colors = c("1" = "red2","2" = "blue3")


###############################################
# Plot of principal component basis functions #
###############################################

ggplot() +
  
  geom_line(mapping = aes(x = eval.grid, y = value, col = factor(PCBF)), data = pcbf.est,   size = 3, show.legend = T, alpha = 1, linetype = 2) +
  geom_line(mapping = aes(x = eval.grid, y = value, col = factor(PCBF)), data = phi_k_list, size = 1.2, show.legend = T, alpha = 1, linetype = 1) +
  
  labs(x = " ", y = " ",title = bquote("n = "* .(n) * ", p = " * .(p.seq[i]) * " and  m = " *.(p.eval+1))) +
  scale_color_manual(values = pcbf.colors, name = "k") +
  
  theme(plot.title   = element_text(size = 22),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.text.x  = element_text(size = 14),     
        axis.text.y  = element_text(size = 14)) +
  
  scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.4, 1.4)) 
ggsave(paste0("Figures/estimation_PCBF.png"), width = 22, height = 14, units = "cm", dpi = 300)




########################################
## Plot of long run covariance kernel ##
########################################

cov.true = cov.matrix(p.eval+1, sigma, theta)

figure1 = plot_ly() |> 
  
  add_markers(x = ~x.eval.grid, y = ~x.eval.grid, z = ~cov.est.1,showlegend = FALSE) |> 
  add_surface(x = ~x.eval.grid, y = ~x.eval.grid, z = ~cov.est.1, showscale = FALSE) |>
  add_surface(x = ~x.eval.grid, y = ~x.eval.grid, z = ~cov.true,  colorscale = list(c(min(cov.true), "#BF382A"),c(max(cov.true), "#0C4B8E")
  ), opacity  = .4, showscale = FALSE) |> 
  
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "", range = c(-0.05, 0.73))))

figure2 = plot_ly() |> 
  
  add_markers(x = ~x.eval.grid, y = ~x.eval.grid, z = ~cov.est.2,showlegend = FALSE) |> 
  add_surface(x = ~x.eval.grid, y = ~x.eval.grid, z = ~cov.est.2, showscale = FALSE) |>
  add_surface(x = ~x.eval.grid, y = ~x.eval.grid, z = ~cov.true,  colorscale = list(c(min(cov.true), "#BF382A"),c(max(cov.true), "#0C4B8E")
  ), opacity  = .4, showscale = FALSE) |> 
  
  layout(scene = list(xaxis = list(title = ""), 
                      yaxis = list(title = ""), 
                      zaxis = list(title = "", range = c(-0.05, 0.73))))

figure1 |> front_layout()
figure2 |> front_layout()


library(reticulate)
use_condaenv("r-reticulate", required = TRUE)
py_config()

save_image(figure1 |> front_layout(), file = "Figures/estimation_covh065.png",   width = 600, height = 750)
save_image(figure2 |> front_layout(),  file = "Figures/estimation_covh025.png",  width = 600, height = 750)
















######################################################################
######################################################################
#####                               ##################################
#####  Optimal bandwidth selection  ##################################
#####                               ##################################
######################################################################
######################################################################



#############################################################
k = 2 # For k = 1,2: functions of analytical solution #######
#############################################################


# Parameters of design.grid and observational rows
p.seq  = c(15,25,50,75,100,150)
n      = 500
p.eval = 200 #(for Simpson rule)

# Setting colors
p.colors = c("15"  = "#F8766D", "25"  = "#A3A500", "50"  = "#00BF7D", "75"  = "#00B0F6", "100" = "#E76BF3", "150" = "purple")



##### optimal bandwidth plot with analytical PCBF ########


#########################
### trapezoidal rule ####
#########################

h.trapez = data.frame()

for(i in 1:6){
  for(kk in 1:k){
    part = paste0("eigfunc_",kk)
    h.trapez = rbind(h.trapez, data.frame(readRDS(paste0("Results/eigenvector/Trapez_p_",p.seq[i],"_n_",n,"_optim.rds"))[[part]]$h_list, k = factor(paste0("k = ",kk)), p = factor(p.seq[i])))
  }
}

ggplot() +
  geom_point(aes(x = h, y = mean.sup, color = p, pch = p), data = h.trapez, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_point(aes(x = h, y = 0.0005 * ( (-1)^seq_along(h) * rep(8, times = length(h)) - 1),color = p, pch = p), data = h.trapez |> group_by(k,p) |> slice_min(mean.sup, n = 1), size = 3,  stroke = 2, alpha = 0.7) +
  scale_color_manual(values = p.colors) +
  labs(x = "h", y = "sup.error",title = bquote("Trapezoidal rule: n = "* .(n) * ", m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/optim_trapezoidal.png"), width = 25, height = 18, units = "cm", dpi = 300)

#####################
### simpson rule ####
#####################

h.simpson = data.frame()

for(i in 1:6){
  for(kk in 1:k){
    part = paste0("eigfunc_",kk)
    h.simpson = rbind(h.simpson, data.frame(readRDS(paste0("Results/eigenvector/Simpson_p_",p.seq[i],"_n_",n,"_optim.rds"))[[part]]$h_list, k = factor(paste0("k = ",kk)), p = factor(p.seq[i])))
  }
}

ggplot() +
  geom_point(aes(x = h, y = mean.sup, color = p, pch = p), data = h.simpson, size = 2.5,  stroke = 1.5, alpha = 0.8, show.legend = F) +
  geom_point(aes(x = h, y = 0.0005 * ( (-1)^seq_along(h) * rep(8, times = length(h)) - 1) ,color = p, pch = p), data = h.simpson |> group_by(k,p) |> slice_min(mean.sup, n = 1), size = 3,  stroke = 2, alpha = 0.7, show.legend = F) +
  
  scale_color_manual(values = p.colors) +
  labs(x = "h", y = "sup.error",title = bquote("Simpson rule: n = "* .(n) * ", m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/optim_simpson.png"), width = 25, height = 18, units = "cm", dpi = 300)



######################################################
### covariance (like in Berger and Holzmann 2025) ####
######################################################

h.cov = data.frame()

for(i in 1:6){
    h.cov = rbind(h.cov, data.frame(readRDS(paste0("Results/covariance/Cov_p_",p.seq[i],"_n_",n,"_optim.rds"))$h_list, k = factor(" "), p = factor(p.seq[i])))
}

ggplot() +
  geom_point(aes(x = h, y = mean.sup, color = p, pch = p), data = h.cov, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_point(aes(x = h, y = 0.0005 * ( (-1)^seq_along(h) * rep(8, times = length(h)) - 1) ,color = p, pch = p), data = h.cov |> group_by(p) |> slice_min(mean.sup, n = 1), size = 3,  stroke = 2, alpha = 0.7) +
  scale_color_manual(values = p.colors) +
  labs(x = "h", y = "sup.error",title = bquote("n = "* .(n) )) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  facet_grid(. ~ k)


ggsave(paste0("Figures/optim_covariance.png"), width = 25, height = 18, units = "cm", dpi = 300)





# eigenvalues


h.trapez = data.frame()

for(i in 1:6){
  for(kk in 1:k){
    part = paste0("eigval_",kk)
    h.trapez = rbind(h.trapez, data.frame(readRDS(paste0("Results/eigenvalue/Trapez_p_",p.seq[i],"_n_",n,"_optim.rds"))[[part]]$h_list, k = factor(paste0("k = ",kk)), p = factor(p.seq[i])))
  }
}

ggplot() +
  geom_point(aes(x = h, y = mean.sup, color = p, pch = p), data = h.trapez, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_point(aes(x = h, y = 0.0005 * ( (-1)^seq_along(h) * rep(8, times = length(h)) - 1),color = p, pch = p), data = h.trapez |> group_by(k,p) |> slice_min(mean.sup, n = 1), size = 3,  stroke = 2, alpha = 0.7) +
  scale_color_manual(values = p.colors) +
  labs(x = "h", y = "sup.error",title = bquote("Trapezoidal rule: n = "* .(n) * ", m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/optim_trapezoidal.png"), width = 25, height = 18, units = "cm", dpi = 300)

#####################
### simpson rule ####
#####################

h.simpson = data.frame()

for(i in 1:6){
  for(kk in 1:k){
    part = paste0("eigval_",kk)
    h.simpson = rbind(h.simpson, data.frame(readRDS(paste0("Results/eigenvalue/Simpson_p_",p.seq[i],"_n_",n,"_optim.rds"))[[part]]$h_list, k = factor(paste0("k = ",kk)), p = factor(p.seq[i])))
  }
}

ggplot() +
  geom_point(aes(x = h, y = mean.sup, color = p, pch = p), data = h.simpson, size = 2.5,  stroke = 1.5, alpha = 0.8, show.legend = F) +
  geom_point(aes(x = h, y = 0.0005 * ( (-1)^seq_along(h) * rep(8, times = length(h)) - 1) ,color = p, pch = p), data = h.simpson |> group_by(k,p) |> slice_min(mean.sup, n = 1), size = 3,  stroke = 2, alpha = 0.7, show.legend = F) +
  
  scale_color_manual(values = p.colors) +
  labs(x = "h", y = "sup.error",title = bquote("Simpson rule: n = "* .(n) * ", m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(-0.01, 0.35)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/optim_simpson.png"), width = 25, height = 18, units = "cm", dpi = 300)

















############################################################################
############################################################################
#####                                     ##################################
#####  supnorm error of for increasing p  ##################################
#####                                     ##################################
############################################################################
############################################################################

#############################################################
k = 2 # For k = 1,2: functions of analytical solution #######
#############################################################

p.seq  = c(15,25,50,75,100,125,150)
n.seq  = c(50,100,200,350,500)
p.eval = 200 #(for Simpson rule)

# Setting colors
n.colors = c("50"  = "#F8766D", "100"  = "#A3A500", "200"  = "#00BF7D", "350"  = "#00B0F6", "500" = "#E76BF3")




######################################
### Eigencector: trapezoidal rule ####
######################################

sup.error.trap = data.frame()

for(i in 1:7){
  for(j in 1:5){
    for(kk in 1:k){
      part = paste0("eigfunc_",kk)
      sup.error      = min(readRDS(paste0("Results/eigenvector/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds"))[[part]]$h_list$mean.sup)
      sup.error.trap = rbind(sup.error.trap, data.frame(sup.error = sup.error, k = factor(paste0("k = ",kk)), p = p.seq[i], n = factor(n.seq[j])))
    }
  }
}




ggplot() +
  geom_point(aes(x = p, y = sup.error, color = n, pch = n),   data = sup.error.trap, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_line( aes(x = p, y = sup.error, color = n, group = n, lty = n), data = sup.error.trap, size = 1) +
  
  scale_color_manual(values = n.colors) +
  labs(x = "p", y = "sup.error",title = bquote("Trapezoidal rule: m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(0, 0.45)) +
  scale_x_continuous(breaks = c(0,50,100,150), limits = c(0,150)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/sup_error_trapezoidal.png"), width = 25, height = 18, units = "cm", dpi = 300)





##################################
### Eigenvector: simpson rule ####
##################################

sup.error.simp = data.frame()

for(i in 1:7){
  for(j in 1:5){
    for(kk in 1:k){
      part = paste0("eigfunc_",kk)
      sup.error      = min(readRDS(paste0("Results/eigenvector/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds"))[[part]]$h_list$mean.sup)
      sup.error.simp = rbind(sup.error.simp, data.frame(sup.error = sup.error, k = factor(paste0("k = ",kk)), p = p.seq[i], n = factor(n.seq[j])))
    }
  }
}


ggplot() +
  geom_point(aes(x = p, y = sup.error, color = n, pch = n),   data = sup.error.simp, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_line( aes(x = p, y = sup.error, color = n, group = n, lty = n), data = sup.error.simp, size = 1) +
  
  scale_color_manual(values = n.colors) +
  labs(x = "p", y = "sup.error",title = bquote("Simpson rule: m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(0, 0.45)) +
  scale_x_continuous(breaks = c(0,50,100,150), limits = c(0,150)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/sup_error_simpson.png"), width = 25, height = 18, units = "cm", dpi = 300)





#####################################
### Eigenvalue: trapezoidal rule ####
#####################################

sup.error.trap = data.frame()

for(i in 1:7){
  for(j in 1:5){
    for(kk in 1:k){
      
      part  = paste0("eigfunc_",kk)
      h.min = as.numeric(readRDS(paste0("Results/eigenvector/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds"))[[part]]$bw.min)
      
      part = paste0("eigval_",kk)
      sup.error      = readRDS(paste0("Results/eigenvalue/Trapez_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds"))[[part]]$h_list |> filter(near(h, h.min)) |> select(mean.sup)
      sup.error.trap = rbind(sup.error.trap, data.frame(sup.error = sup.error$mean.sup, k = factor(paste0("k = ",kk)), p = p.seq[i], n = factor(n.seq[j])))
    }
  }
}


ggplot() +
  geom_point(aes(x = p, y = sup.error, color = n, pch = n),   data = sup.error.trap, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_line( aes(x = p, y = sup.error, color = n, group = n, lty = n), data = sup.error.trap, size = 1) +
  
  scale_color_manual(values = n.colors) +
  labs(x = "p", y = "error",title = bquote("Trapezoidal rule: m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(0, 0.45)) +
  scale_x_continuous(breaks = c(0,50,100,150), limits = c(0,150)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/error_EV_trapezoidal.png"), width = 25, height = 18, units = "cm", dpi = 300)




#################################
### Eigenvalue: simpson rule ####
#################################

sup.error.simp = data.frame()

for(i in 1:7){
  for(j in 1:5){
    for(kk in 1:k){
      
      part  = paste0("eigfunc_",kk)
      h.min = as.numeric(readRDS(paste0("Results/eigenvector/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds"))[[part]]$bw.min)
      
      part = paste0("eigval_",kk)
      sup.error      = readRDS(paste0("Results/eigenvalue/Simpson_p_",p.seq[i],"_n_",n.seq[j],"_optim.rds"))[[part]]$h_list |> filter(near(h, h.min)) |> select(mean.sup)
      sup.error.simp = rbind(sup.error.simp, data.frame(sup.error = sup.error$mean.sup, k = factor(paste0("k = ",kk)), p = p.seq[i], n = factor(n.seq[j])))
    }
  }
}


ggplot() +
  geom_point(aes(x = p, y = sup.error, color = n, pch = n),   data = sup.error.simp, size = 2.5,  stroke = 1.5, alpha = 0.8) +
  geom_line( aes(x = p, y = sup.error, color = n, group = n, lty = n), data = sup.error.simp, size = 1) +
  
  scale_color_manual(values = n.colors) +
  labs(x = "p", y = "error",title = bquote("Simpson rule: m = " * .(p.eval+1))) +
  
  theme(plot.title   = element_text(size = 25),
        legend.text  = element_text(size = 22),
        legend.title = element_text(size = 22),
        strip.text   = element_text(size = 22),
        axis.title.x = element_text(size = 22),     
        axis.title.y = element_text(size = 22),
        axis.text.x  = element_text(size = 20),     
        axis.text.y  = element_text(size = 20),
        
        legend.key.size = unit(0.8, "cm")) +
  
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(limits = c(0, 0.45)) +
  scale_x_continuous(breaks = c(0,50,100,150), limits = c(0,150)) +
  facet_grid(. ~ k)

ggsave(paste0("Figures/error_EV_simpson.png"), width = 25, height = 18, units = "cm", dpi = 300)
