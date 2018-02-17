library(R.matlab)
library(svd)
library(corpcor) # covariance matrix estimation package
library(ggplot2)
library(Kendall)
library(pracma)

rm(list = ls())
# I obtained the following data from my Matlab simulation
# It is an n by m matrix where n is the number of states and m is the number of snapshots
# raw_data <- readMat('harvesting_moving_new.mat') 
raw_data <- read.csv('veg_turb.csv')
# raw_data <-  do.call(rbind, raw_data)
norm_vec <- function(x) sqrt(sum(x^2))

# define hyperparameters such as window size, moving step
inter <- 6
wind_size <- 3000
# moving_win <- seq(from = 1, to = 1000, by = 1000)
moving_win <- seq(from = 1, to = 40000-wind_size, by = 1000)

svd_c_sample <- rep(1,length(moving_win))
svd_c_shrink <- rep(1,length(moving_win))
svd_c_shrink_per <- rep(1,length(moving_win))
k <- 0

# perform covariance matrix estimation using sample covariance matrix estimation and 
# the shrinkage approach
for (i in moving_win){
k <- k + 1  
data <- t(raw_data[seq(from = i, to = i+wind_size, by = inter),])

# detrend
# M <- rowMeans(data,na.rm = FALSE, dims = 1)
# M <- matrix(rep(M,times=dim(data)[2]),nrow=dim(data)[1],ncol=dim(data)[2])
# data <- data - M
data <- detrend(data, tt='linear')
# sample covariance matrix estimation
s1 <- cov(t(data))
s1_svd <- svd(s1)
svd_c_sample[k] <- s1_svd[[1]][[1]]
# shrinkage covariance matrix estimation
s2 <- cov.shrink(t(data))
s2_svd <- svd(s2)
svd_c_shrink[k] <- s2_svd[[1]][[1]]
svd_c_shrink_per[k] <- s2_svd[[1]][[1]]/norm_vec(s2_svd[[1]])
}

write.csv(svd_c_shrink, 'veg_turb_shrink.csv', row.names = FALSE)
write.csv(svd_c_shrink_per, 'veg_turb_shrink_per.csv', row.names = FALSE)