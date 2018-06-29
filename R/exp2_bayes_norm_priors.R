###############################################################################
#																			  #
# execute exp2_bayes_norm_priors.R									          #
#																			  #
###############################################################################

setwd("Documents/wiskunde/2017-2018/bachelor_project/R/handin")
setwd("/media/mynewdrive1/Documenten/Wiskunde/2017-2018/bachelor_project/R/handin")

# libraries
library("arm")
# load functions
source("helpers.R")
source("helpers_CV_bayesglm.R")
# load variables
source("variables_exp2.R")

# Choose what initial betas and sample size to use
## Betas
k = 1
## Sample size
l = 1

# for testing
repeat_cross_val <- 1001
alpha <- seq(0.001,15,length=101)
order_fit <- 2
alpha_plot =  TRUE; data_plot = TRUE
alpha_plot =  FALSE; data_plot = FALSE
orthog <- TRUE
orthog <- FALSE

# Execute experiment
alpha_stars <- c()
mses <- c()
mu_news <- c()
for (i in 1:repeat_cross_val){
  cross_val_temp <- execute_cross_val(alpha = alpha, x_min = x_min, x_max = x_max, order_sample = order_p, order_fit = order_fit, initial_betas = initial_betas[[k]], sample_size = sample_size[l], orthog =orthog, alpha_plot =  alpha_plot, data_plot = data_plot, multiple_samples = FALSE)
  #
  alpha_stars[i] <- cross_val_temp[[1]]$alpha_star
  mses[i] <- cross_val_temp[[1]]$alpha_star_cross_val$mse
  mu_news[[i]] <- cross_val_temp[[1]]$alpha_star_cross_val$mu_new
}

# calc out of back mse = test mse using mu_news
i <- 1
test_mses <- c()
for (beta in mu_news[1:1001]) {
  test_mses[i] <- out_of_back_mse(beta_hat = beta, beta_true = c(1,1,1), x_min = x_min, x_max = x_max, n = 1000)
  i <- i + 1
}
mean(test_mses)

# Visualise experiment
## alphas
hist(alpha_stars,xlab=paste(expression(alpha),"*",sep=''),ylab="Frequentie",main="")
plot(alpha_stars,mses,xlab=paste(expression(alpha),"*",sep=''),ylim=c(0,13))
plot(alpha_stars)
## (test) mse's
hist(mses,breaks =20,main="",xlab="MSE",ylab="Frequentie")
hist(test_mses,main="",xlab="MSE",ylab="Frequentie")
hist(sort(test_mses)[1:(0.9*length(test_mses))],breaks=20,main="",xlab="MSE",ylab="Frequentie")


# numbers
mean(alpha_stars)
mean(mses)
mu_1 <- 0
mu_2 <- 0
mu_3 <- 0
mu_4 <- 0
mu_5 <- 0
mu_6 <- 0
for (j in 1:length(mu_news)){
  mu_1 <- mu_1 + mu_news[[j]][1]
  mu_2 <- mu_2 + mu_news[[j]][2]
  mu_3 <- mu_3 + mu_news[[j]][3]
  mu_4 <- mu_4 + mu_news[[j]][4]
  mu_5 <- mu_5 + mu_news[[j]][5]
  mu_6 <- mu_6 + mu_news[[j]][6]
}
mu_1_mean <- mu_1 / length(mu_news)
mu_2_mean <- mu_2 / length(mu_news)
mu_3_mean <- mu_3 / length(mu_news)
mu_4_mean <- mu_4 / length(mu_news)
mu_5_mean <- mu_5 / length(mu_news)
mu_6_mean <- mu_6 / length(mu_news)

mu_1_mean
mu_2_mean
mu_3_mean
mu_4_mean
mu_5_mean
mu_6_mean


