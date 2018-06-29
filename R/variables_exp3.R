###############################################################################
#																			  #
# variables for exp3_bayes_tpriors.R									      #
#																			  #
###############################################################################

# polynomial order of x from which sample y = b_0 + b_1 x + b_2 x^2 + ... + b_order_p x^order_p is drawn
order_p <- 2
#polynomial order of x which is assumed in the model (equal or different from order_p)
order_fit <- 2
# values of x in [x_min, x_max]
x_min <- -1
x_max <- 3
# sample sizes
sample_size <- c(20,60,300)

# Cross validation
# number_of_partitions = folds
number_of_partitions <- 5 
# number of times cross-validations has to be executed with same data before moving on to new data (with new conditions)
repeat_cross_val <- 1001

# save all mses
mse_list <- list()
# save mu_news
mu_news <- list()

# initial betas
initial_betas <- list(seq(1,1,length=order_fit + 1),seq(10,10,length=order_fit + 1),seq(1,1,length=order_fit + 1),c(1,0.1,0.01,0.001,0.0001,0.00001))

# different degrees of freedom
df <- seq(1,100,length=100)
df <- seq(1,100,by=3)

#
df_stars <- c()
mses <- c()
