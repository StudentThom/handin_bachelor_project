###############################################################################
#																			  #
# execute exp3_bayes_t_priors.R									              #
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
source("variables_exp3.R")

# Choose what initial betas and sample size to use
## Betas
k = 1
initial_betas[[k]] = c(1,1,1,1,1,1)
## Sample size
l = 1
order_fit <- 2
repeat_cross_val <-1001
df_plot =  TRUE; data_plot = TRUE
df_plot =  FALSE; data_plot = FALSE

# Execute experiment
for (i in 1:repeat_cross_val){
  cross_val_temp <- execute_cross_val_tprior(df = df, x_min = x_min, x_max = x_max, order_sample = order_p, order_fit = order_fit, initial_betas = initial_betas[[k]], sample_size = sample_size[l], orthog =FALSE, df_plot =  df_plot, data_plot = data_plot, multiple_samples = FALSE)
  #
  df_stars[i] <- cross_val_temp[[1]]$df_star
  mses[i] <- cross_val_temp[[1]]$df_star_cross_val$mse
  mu_news[[i]] <- cross_val_temp[[1]]$df_star_cross_val$mu_new
  print(i)
}

# calc out of back mse = test mse using mu_news
i <- 1
test_mses <- c()
for (beta in mu_news[1:1001]) {
  #print(beta)
  test_mses[i] <- out_of_back_mse(beta_hat = beta, beta_true = c(1,1,1), x_min = x_min, x_max = x_max, n = 100)
  i <- i + 1
}
mean(test_mses)
mean(test_mses,trim=0.1)
hist(test_mses)
hist(test_mses,breaks=10000,main="",xlab="MSE",ylab="Frequentie")
hist(sort(test_mses)[1:(0.8*length(test_mses))],breaks=10,main="",xlab="MSE",ylab="Frequentie")

# Visualise experiment
hist(df_stars,main="",xlab="df",ylab="Frequentie")
plot(df_stars,mses)
plot(df_stars)
hist(mses,breaks =50,main="",xlab="MSE",ylab="Frequentie")
# numbers
mean(df_stars)
mean(mses)
mean(mses,trim=0.1)
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

mu_news_mean <- c(mu_1_mean,mu_2_mean,mu_3_mean,mu_4_mean,mu_5_mean,mu_6_mean)
# bias?!
verschil1 <- abs(c(1,1,1,1,1,1)-c(mu_1_mean,mu_2_mean,mu_3_mean,mu_4_mean,mu_5_mean,mu_6_mean))
bias_sq <- mean(verschil1)**2
bias_sq

verschil2 <- 0
var1 <- 0
for (beta in mu_news[1:1000]){
  # bias
  #verschil2 <- verschil2 + abs(c(beta[1,1],beta[2,1],beta[3,1]) - c(1,1,1))
  verschil2 <- verschil2 + abs(c(beta[1],beta[2],beta[3]) - c(1,1,1))
  #print(verschil2)
  # var
  #var1 <- var1 + sum((mu_news_mean - beta)**2)
}
bias_sq2 <- mean(verschil2 / 1000)**2
bias_sq2
var1 <- var1 / 1000
var1

# Var?
var1 <- mean(mean(mu_news)-mu_news)

