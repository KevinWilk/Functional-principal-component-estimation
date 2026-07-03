###############################################################################################################################################################################################################################################################
######  Estimates from Github repository github.com/KevinWilk/Beyond-average-warming ##########################################################################################################################################################################
###############################################################################################################################################################################################################################################################

# Loading estimated lag covariance kernels for lag = 0,...,12 from Wilk and Holzmann (2026)
tmp          = tempfile(fileext = ".rds")
download.file("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/weather%20temperature/long%20run%20kernel/Results/list_Gamma_d_Berlin.rds",
              destfile = tmp, mode = "wb")
list.Gamma.d = readRDS(tmp)

# Loading maximum lag of long run kernel from Wilk and Holzmann (2026)
tmp          = tempfile(fileext = ".rds")
download.file("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/weather%20temperature/long%20run%20kernel/Results/max_lag_Berlin.rds",
              destfile = tmp, mode = "wb")
max.lag      = readRDS(tmp)

rm(tmp)



# Calculating long run covariance matrix for each month from Wilk and Holzmann (2026)

start.24 = 21; end.24 = 117 

lr.Gamma.d.month = lapply(1:12, function(m) {lr.cov  = Reduce(`+`,lapply(1:(max.lag[m,2]+1),  
                                                                         function(j){if(j >= 2){
                                                                           (list.Gamma.d[[m]][[j]][start.24:end.24,start.24:end.24]+t(list.Gamma.d[[m]][[j]][start.24:end.24,start.24:end.24]))*(1-(j-1)/(max.lag[m,2]+1))
                                                                         }else{list.Gamma.d[[m]][[j]][start.24:end.24,start.24:end.24]} 
                                                                         }))
return(lr.cov)})


# Loading data set of Berlin 2000 - 2025 from Wilk and Holzmann (2026)
tmp          = tempfile(fileext = ".RData")
download.file("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/weather%20temperature/data%20sets/Berlin.RData",
              destfile = tmp, mode = "wb")
load(tmp)

# Loading weight of bivariate kernel estimator for lag = 0 from Wilk and Holzmann (2026)
tmp          = tempfile(fileext = ".rds")
download.file("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/weather%20temperature/kernel%20weights/full_w_d_lag0_006.rds",
              destfile = tmp, mode = "wb")
w0 = readRDS(tmp)


### January of Berlin 2000 - 2025
jan.data  = data.d.34h |> 
  filter(MONTH == month.name[1]) |> 
  dplyr::select(4:dim(data.d.34h)[2])

jan.data = jan.data[rowSums(is.na(jan.data)) == 0,]

### estimated lagged covariance from lag 0 to 2
jan.cov.list = list.Gamma.d[[1]][1:(max.lag[1,2]+1)]
jan.cov.lag0 = jan.cov.list[[1]][start.24:end.24,start.24:end.24]

### August of Berlin 2000 - 2025
aug.data  = data.d.34h |> 
  filter(MONTH == month.name[8]) |> 
  dplyr::select(4:dim(data.d.34h)[2])

aug.data = aug.data[rowSums(is.na(aug.data)) == 0,]

### estimated lagged covariance from lag 0 to 2
aug.cov.list = list.Gamma.d[[8]][1:(max.lag[8,2]+1)]
aug.cov.lag0 = aug.cov.list[[1]][start.24:end.24,start.24:end.24]

# remove data sets from Wilk and Holzmann (2026)
rm(data.d.24h.hourly,data.d.30h, data.d.34h, data.s.24h, data.s.34h, list.Gamma.d)
rm(end.24, start.24, filename, tmp)

save.image("weather temperature/Berlin_lag0_setup.RData")



###############################################################################################################################################################################################################################################################



load("weather temperature/Berlin_lag0_setup.RData")



# Loading weight of bivariate kernel estimator for lag = 1,2 from Wilk and Holzmann (2026)

tmp          = tempfile(fileext = ".rds")
download.file("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/weather%20temperature/kernel%20weights/full_w_d_lag1_006.rds",
              destfile = tmp, mode = "wb")
w1 = readRDS(tmp)

tmp          = tempfile(fileext = ".rds")
download.file("https://raw.githubusercontent.com/KevinWilk/Beyond-average-warming/main/weather%20temperature/kernel%20weights/full_w_d_lag2_006.rds",
              destfile = tmp, mode = "wb")
w2 = readRDS(tmp)

rm(tmp,filename)

save.image("weather temperature/Berlin_lagb_setup.RData")

