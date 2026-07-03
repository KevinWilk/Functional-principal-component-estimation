######################################################################
####################### function  ####################################
######################################################################

mu0 = function(x){
  return(rep(0, length(x)))
}

#################################################
################  Ornstein-Uhlenbeck  ###########
#################################################


OU = function(n, t = seq(0, 1, len = 201), mu = 0, alpha = 1, sigma = 1, stat = T){
  
  if(stat){x0 = rnorm(n, mean = 0,sd = sigma/sqrt(2 * alpha))}
  else{    x0 = 0}
  
  goffda::r_ou(n, t, mu, alpha, sigma, x0)$data
}

cov_ou = function(t, sigma, theta){
  sigma^2/(2*theta) *
    exp(-theta * abs(t[1] - t[2]))
}


cov.matrix = function(p.eval, sigma = 1, theta = 1){
  
  M = matrix(0, p.eval, p.eval)
  M.up = upper.tri(M, T)
  
  obs.grid = function(p = NULL, x = NULL, comp = "less") {
    
    if(is.null(x) & is.null(p))
      stop("p or grid needs to be supplied")
    else if(is.null(x))
      x = 0:(p-1)/(p-1)
    x.grid = expand.grid(x, x)
    
    if (comp == "less")
      b = x.grid[,1] < x.grid[,2]
    if (comp == "lesseq")
      b = x.grid[,1] <= x.grid[,2]
    if (comp == "gtr")
      b = x.grid[,1] > x.grid[,2]
    if (comp == "gtreq")
      b = x.grid[,1] >= x.grid[,2]
    if (comp == "without diagonal")
      b = !(x.grid[,1] == x.grid[,2])
    if (comp == "full")
      b = T
    x.grid[b, ]
  }
  
  M[M.up] = cov_ou(obs.grid(p.eval, comp = "lesseq"), sigma, theta)$Var1 
  M[!M.up] = t(M)[!M.up]
  return(M)
}

#################################################
### functions for analytical solution ###########
#################################################

# solve euqation (...) to otain omega_k
cond = function(omega,theta = 1){
  (theta^2-omega^2) * sin(omega) + 2 * theta * omega * cos(omega)
}

omega_k = function(k, theta){
  uniroot(
    function(omega) cond(omega, theta),
    interval = c((k - 1) * pi + 1e-8, k * pi - 1e-8)
  )$root
}

# analytical principal component basis function
phi_k = function(x, omega = 0, theta = 1){
  
  Ak = (integrate(function(y, w = omega, a = theta){( cos(w * y) + a/w * sin(w * y) )^2 },lower = 0, upper = 1)$value)^{-1/2}
  
  return( Ak * (cos(omega * x) + theta/omega * sin(omega * x))  )
}


################################################
######## Qadrature formulas ####################
################################################

# Composite trapezoidal rule
trapezoid_weights = function(p.eval){
  1/p.eval * c(1/2,rep(1,times = p.eval -2),1/2) 
}

# Composite Simpson rule
simpson_weights = function(p.eval){
  p = p.eval - 1
  if(p %% 2 != 0){stop("p.eval-1 must be even")}
  w = rep(c(4,2), length.out = p - 1)
  1/(3*p.eval) * c(1, w, 1)
}




################################################################################################
########### bandwidth selection via K-fold cross validation ####################################
################################################################################################



h.optim = function(N, n, p, p.eval, h.seq, 
                   k = 1, m = 1, w.parallel = T, 
                   f, r.process, process.arg, eps.arg,
                   quadrature = "trapezoidal", supnorm = "eigenfunction"){
  
  if(quadrature %in% c("trapezoidal","Simpson")){
    
    x.design    = 0:p/p
    x.eval.grid = 0:p.eval/p.eval
    
    help = function(h){
      w_h      = local.polynomial.weights(p + 1, h, p.eval = p.eval + 1, m = m, parallel = w.parallel, parallel.environment = F, grid.type = "less") # from Wilk and Holzmann (2026)
      max_diff = numeric(1)
      
      omega_list  = sapply(1:k, omega_k, theta = theta)
      lambda_list = sapply(omega_list, function(x, sigma, theta){sigma^2/(theta^2+x^2)}, sigma = sigma, theta = theta)
      phi_k       = sapply(omega_list, function(x, omega, theta){ phi_k(x , omega, theta)  }, theta = theta, x = x.eval.grid)[,k]
      
      future_replicate(N, {
        Y        = FDA_observation(n, x.design, f, r.process = r.process, process.arg = process.arg, eps.arg = eps.arg)
        cov.est  = eval.weights(w_h, observation.transformation(Y, grid.type = "less"))
        #cov.true = cov.matrix(p.eval+1, sigma, theta)
        
        # Quadrature formula
        if( quadrature == "trapezoidal"){     quad.w = trapezoid_weights(p.eval+1) }
        else if(quadrature == "Simpson"){     quad.w = simpson_weights(p.eval+1)   }
        
        # calulate pcb function
        result  = eigen( diag(sqrt(quad.w)) %*% cov.est %*% diag(sqrt(quad.w)), symmetric = TRUE) 
        eig.val = result$values[k]
        
        #result.true  = eigen( diag(sqrt(quad.w)) %*% cov.true %*% diag(sqrt(quad.w)), symmetric = TRUE) 
        #eig.val.true = result.true$values[k]
        #eig.vec.true = result.true$vectors[,k] * sqrt(quad.w)
        
        sign    = as.vector(sign(crossprod(result$vectors[,k], phi_k)))
        eig.vec = result$vectors[,k] * sqrt(quad.w)
        
        eig.vec.est = cov.est   %*% ( eig.vec  * 1/eig.val * sign) 
        
        if(     supnorm == "eigenvalue"){    max_diff = abs((eig.val - lambda_list[k])/lambda_list[k]) }
        else if(supnorm == "eigenfunction"){ max_diff = max(abs(eig.vec.est - phi_k) )}
        else{stop("supnorm difference not implemented")}
        
        max_diff
      }, future.seed = T)
    }
    
    mean_sup = sapply(h.seq, help)
    
    h_list = data.frame(h = h.seq, mean.sup = apply(mean_sup, 2, mean))
    
    return(list(bw.min = h.seq[which.min(apply(mean_sup, 2, mean))], h_list = h_list)) 
    
  }else{stop("Quadrature formula not implemented")}
}



h.optim.cov = function(N, n, p, p.eval, h.seq, 
                       m = 1, w.parallel = T, 
                       f, r.process, process.arg, eps.arg){ 
  
  x.design    = 0:p/p
  cov.true = cov.matrix(p.eval+1, sigma, theta)  
  
  help = function(h){
    w_h      = local.polynomial.weights(p + 1, h, p.eval = p.eval + 1, m = m, parallel = w.parallel, parallel.environment = F, grid.type = "less") # from Wilk and Holzmann (2026)
    max_diff = numeric(1)
    
    future_replicate(N, {
      Y        = FDA_observation(n, x.design, f, r.process = r.process, process.arg = process.arg, eps.arg = eps.arg)
      cov.est  = eval.weights(w_h, observation.transformation(Y, grid.type = "less"))
      
      max_diff = max(abs(cov.true - cov.est) )
      
      max_diff
    }, future.seed = T)
  }
  
  mean_sup = sapply(h.seq, help)
  
  h_list = data.frame(h = h.seq, mean.sup = apply(mean_sup, 2, mean))
  
  return(list(bw.min = h.seq[which.min(apply(mean_sup, 2, mean))], h_list = h_list)) 
}



### Multiplier Bootstrap


mb.PCBF = function(N, sample, list.weights, cov.est, quadrature = "trapezoidal", k.max = 1, dependent = F, periodic = F, m = NA, max.lag = 0){
  
  if(length(list.weights) != (max.lag + 1) ){ stop(paste0("length of long-run kernel weights list must be ",max.lag + 1)) }
  
  n      = dim(sample)[1]
  p      = dim(sample)[2]
  p.eval = dim(cov.est)[2]
  
  if(quadrature %in% c("trapezoidal","Simpson")){
    
    # Quadrature formula
    if( quadrature == "trapezoidal"){     quad.w = trapezoid_weights(p.eval) }
    else if(quadrature == "Simpson"){     quad.w = simpson_weights(p.eval)   }
    
    #calculating PCBF and eigenvalue
    result.est  = eigen( diag(sqrt(quad.w)) %*% cov.est  %*% diag(sqrt(quad.w)), symmetric = T) 
    eig.val.est = result.est$values[1:k.max]
    eig.vec.est = result.est$vectors[,1:k.max] * sqrt(quad.w) 
    eig.vec.est = as.matrix(cov.est   %*% ( eig.vec.est  %*% diag(1/eig.val.est)) ) 
    
    #multiplier bootstrap repetition
    sample.mb = future_replicate(N, {
      
      if(dependent == T){
        
        
        ln_func = function(n){floor(2*n^(1/3))}
        k1 = function(h, n, func) {
          L = func(n)
          ifelse(abs(h) < L, 1 / (2*L - 1), 0)
        }
        
        
        q_n = 1/(2*ln_func(n)-1)
        w_n  = rnorm(3*n, mean = 0, sd = 1/sqrt(q_n))
        g_n  = numeric(n)
        
        for(j in 1:n){
          g_n[j] = sum(sapply((-ln_func(n)):ln_func(n), function(h) k1(h, n, ln_func)) * w_n[j:(j+2*ln_func(n))])
        }
        
        g_n  = g_n  - mean(g_n)
        
      }else{g_n  = rnorm(n)}
      
      
      # Multiplier Bootstrap bei cov.est 
      
      if(periodic == T){
        
        if(is.na(m)){m = p - 1}
        map = ((1:p) - 1) %% m + 1
        
        sample.means.base = vapply(1:m,function(g) mean(unlist(sample[, which(map == g)]), na.rm = T),numeric(1))
        sample.means = sample.means.base[map]
        
      }else{
        sample.means = colMeans(sample, na.rm = T)  
      }
      
      # lag 0 kernel terms #####################################################################################################################################
      sample.2    = do.call(cbind, lapply( lapply(1:n, function(i) tcrossprod(unlist( sample[i,]), unlist(sample[i,])) ), as.vector)) 
      multipliers = do.call(cbind,lapply(1:n, function(i)as.vector(tcrossprod(g_n[i] * (numeric(p) + 1), numeric(p) + 1))))                                         
      Z           = (sample.2 - as.vector(tcrossprod(sample.means))) * multipliers   
      
      Z.mean =  Z |>
        apply(1, function(x){
          sum(x, na.rm = T)/(n - 1 - sum(is.na(x)))})  
      
      if(periodic == T){
        
        col = rep(1:p, each = p)
        row = rep(1:p, times = p)
        group = (map[row] - 1) * m + map[col]
        Z.mean.base = tapply(Z.mean, group, mean, na.rm = T)
        Z.mean = as.numeric(Z.mean.base[group])
        
      }
      
      # Using "less" as grid for estimating kernels lag 0
      M      = as.vector(upper.tri(matrix(0, p, p)))  
      cov.mb = eval.weights(list.weights[[1]], Z.mean[M])
      start.24 = 21; end.24 = 117 # Not optimal implemented
      cov.mb = cov.mb[start.24:end.24,start.24:end.24] # reduce to 24hx24h Matrix 
      if( any(dim(cov.mb) != dim(cov.est)) ){stop("Dimensions of estimated kernel und kernel with multipliers have to be equal")}                                         
      
      # lag b>0 kernel term ####################################################################################################################################
      if(max.lag > 0){
        for(b in 1:max.lag){
          sample1.lag = do.call(cbind, lapply(lapply(1:(n - abs(b)), function(i) tcrossprod(unlist(sample[i,]), unlist(sample[i+abs(b),]))), as.vector))     
          sample2.lag = do.call(cbind, lapply(lapply((abs(b)+1):n,   function(i) tcrossprod(unlist(sample[i,]), unlist(sample[i-abs(b),]))), as.vector))
          
          multipliers.1 = do.call(cbind,lapply(1:(n - abs(b)), function(i)as.vector(tcrossprod(g_n[i] * (numeric(p) + 1), numeric(p) + 1))))
          multipliers.2 = do.call(cbind,lapply((abs(b)+1):n,   function(i)as.vector(tcrossprod(g_n[i] * (numeric(p) + 1), numeric(p) + 1))))
          
          Z1 = ( sample1.lag - as.vector(tcrossprod(sample.means)) ) * multipliers.1      
          Z2 = ( sample2.lag - as.vector(tcrossprod(sample.means)) ) * multipliers.2                                           
          
          Z1.mean =  Z1 |>
            apply(1, function(x){
              sum(x, na.rm = T)/(n - 1 - sum(is.na(x)))})                                        
          Z2.mean =  Z2 |>
            apply(1, function(x){
              sum(x, na.rm = T)/(n - 1 - sum(is.na(x)))})
          
          if(periodic == T){
            col = rep(1:p, each = p)
            row = rep(1:p, times = p)
            group = (map[row] - 1) * m + map[col]
            
            Z1.mean.base = tapply(Z1.mean, group, mean, na.rm = T)
            Z1.mean =  as.numeric(Z1.mean.base[group])
            
            Z2.mean.base = tapply(Z2.mean, group, mean, na.rm = T)
            Z2.mean =  as.numeric(Z2.mean.base[group])
          }
          
          # Using "lesseq" as grid for estimating kernels lag b > 0                                               
          M.lag      = as.vector(upper.tri(matrix(0, p, p)) | diag(T, p, p))
          cov.lag.mb = eval.weights(list.weights[[b+1]], list(Z1 = Z1.mean[M.lag], Z2 = Z2.mean[M.lag], lag = b), lag = b)
          cov.lag.mb  = cov.lag.mb[start.24:end.24,start.24:end.24] # reduce to 24hx24h Matrix  
          
          # Long run kernel estimator from Wilk and Holzmann (2026): weights from Newey and West (1994) 
          cov.mb = cov.mb + (1 - (b/(max.lag + 1)) ) * ( cov.lag.mb + t(cov.lag.mb) )                                                
        }
      }
      
      #calculating PCBF and eigenvalue
      result.mb  = eigen( diag(sqrt(quad.w)) %*% cov.mb  %*% diag(sqrt(quad.w)), symmetric = T) 
      eig.val.mb = result.mb$values[1:k.max]
      eig.vec.mb = result.mb$vectors[,1:k.max] * sqrt(quad.w)
      sgn        = sign(diag(crossprod(eig.vec.mb, eig.vec.est)))
      eig.vec.mb = cov.mb  %*%  ( eig.vec.mb  %*% diag(1/eig.val.mb) %*% diag(sgn) )
      
      
      list(eig.value = eig.val.mb, PCBF = eig.vec.mb / sqrt(colSums(eig.vec.mb^2 * trapezoid_weights(p.eval) )) )
      
    }, future.seed = T)
    
    mb.EV           = as.data.frame(do.call(rbind, sample.mb[1, ]))
    colnames(mb.EV) = paste0("EV", 1:ncol(mb.EV))
    
    mb.vec = as.list(as.data.frame(do.call(rbind, sample.mb[2, ])))
    mb.vec = lapply(mb.vec, function(x) {as.data.frame(matrix(x, nrow = p.eval, byrow = F))})
    
    return( list(eig.value = mb.EV, eig.vec = mb.vec) )
    
  }else{stop("Quadrature formula not implemented")}
}

