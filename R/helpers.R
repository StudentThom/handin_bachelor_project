###############################################################################
#' helpers.R 																                                  #
#' helpers file for basic and polynomial bayesian regression                  #
#' Created by StudentThom, spring 2018										                    #
###############################################################################
#'
#' In this file:
#'
#' create_matrix_x(x, order_p)

#' generate_polynomial_sample(order_p =2, slopes=c(1,1,1), size=100, mean=0, sd=0.2, x_min=0,x_max=1, plaatjes=FALSE)

#' ordinair_least_squres_poly(x, y, order_p = 2, x_min = 0, x_max = 1, plaatjes=FALSE)

#' overfit(order_sample, order_fit, slopes = NULL ,size=50, x_min = 0, x_max = 2, number_of_test_data_samples = 1, plaatjes=TRUE, x_range = NULL, y_range = NULL, experiment_number = 99)

#' plot_sample(x, y, beta_hat, slopes, x_range, y_range, sample_line = FALSE)

#' cross_validation(y, matrix_x, number_of_partitions, alpha)

#' calculate_mse(y, y_hat){

#' make_training_and_test_data_x(number_of_partitions, partitions, leave_out)

#' make_training_and_test_data_y(number_of_partitions, partitions, leave_out)

#' posterior_normal(matrix_x, y, alpha, beta_error)

#' execute_cross_val(alpha, x_min = -1, x_max = 1, order_sample = 2, order_fit = 2, initial_betas=list(c(1,1,1)), sample_size = c(20), orthog = FALSE, alpha_plot = TRUE, data_plot = TRUE, multiple_samples = FALSE)

#' out_of_back_mse(beta_hat, x_min, x_max, n = 100)

#'
#' Make matrix of a vector x. On the rows the entries of x. 
#' On the columns the different orders of x:
#' First column intercept (==1), second column x, third column x^2, ..., x^order_p
#'
create_matrix_x <- function(x, order_p){
  if (order_p <= 0){
    return("Please give order_p higher then zero")
  }
  size <- length(x)
  matrix_x <- matrix(0,nrow=size,ncol=order_p+1)
  matrix_x[,1] = 1
  for (j in 1:order_p){
    matrix_x[,j+1] <- x^j
  }
  return(matrix_x)
}

#'
#' Generate random sample of x ~ unif(xmin, xmax, n)
#' Generate Y ~ B_0 + B_1 X + ... + B_{order_p} x^{order_p} + epsilon, epsilon ~ Norm(mu,sigma^2)
#' Return list with x as matrix (in form of create_matrix_x()) and y as vector
#' 
generate_polynomial_sample <- function(order_p =2, slopes=c(1,1,1), size=100, mean=0, sd=0.2, x_min=0,x_max=1, plaatjes=FALSE){
  # prevent errors
  if (length(slopes) <= order_p){
    print("Please provide enough slopes for the order you use")
    return("Please provide enough slopes for the order you use")
  }
  # Generate sample y = slope * x, x ~ unif[x_min,x_max]
  x <- runif(x_min,x_max,n=size)
  y <- seq(0,0,length=size)
  matrix_x <- create_matrix_x(x, order_p)
  for (i in 1:length(y)){
    for (j in 1:(order_p+1)){
      #error <- rnorm(n=1,mean=mean,sd=sd)
      y[i] <- (y[i] + (slopes[j]*(matrix_x[i,j]))) #+ error
    }
    error <- rnorm(n=1,mean=mean,sd=sd)
    y[i] <- y[i] + error
  }
  if (plaatjes){
    ## plaatjes
    #hist(y)
    plot(x,y,ylim=c(min(y)-0.1,max(y)+0.1))
  }
  return(list("matrix_x"=matrix_x,"y"=y))
}

#'
#' Given a dataset (X,Y), give OLS estimator B for Y = B %*% X
#' Return B == beta_hat
#' 
ordinair_least_squres_poly <- function(x, y, order_p = 2, x_min = 0, x_max = 1, plaatjes=FALSE){
  matrix_x <- create_matrix_x(x,order_p = order_p)
  transposed <- t(matrix_x)
  beta_hat <- solve(transposed%*%matrix_x) %*% transposed %*%y
  
  #plot
  if (plaatjes){
    y_hat <- seq(0,0,length=dim(matrix_x)[1])
    x_seq <- seq(x_min,x_max,length=dim(matrix_x)[1])
    x_seq_matrix <- create_matrix_x(x_seq, order_p = length(beta_hat))
    for (n in 1:length(y_hat)){
      for (p in 1:length(beta_hat))
        #y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq[n]^p
        y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq_matrix[n,p]
    }
    plot(matrix_x[,2],y,ylim=c(0,max(y)+0.1),xlab="x",ylab="y")
    lines(x_seq,y_hat,col="orange")
  }
  return(beta_hat)
}

#'
#' Give an example of overfitting using OLS estimator.
#' An sample of order_sample is draw. Then it is fitted using OLS of order order_fit 
#' (usualy order_sample < order_fit) to demonstrate overfitting.
#' Then the 'overfitted' estimates are used to model a new sample of order_sample, 
#' which will have much higher MSE. This proves the point of overfitting.
#' Ggrahpics for extra explanations are optional.
#'
overfit <- function(order_sample, order_fit, slopes = NULL ,size=50, x_min = 0, x_max = 2, number_of_test_data_samples = 1, plaatjes=TRUE, x_range = NULL, y_range = NULL, experiment_number = 99){
  # write output to new directory
  # if (dir.exists(paste("plaatjes/poly_regr/overfit/exp_",experiment_number,sep = ""))){
  #   print("about to overwrite an older experiment")
  #   return("error, experiment number already exists")
  # }
  # else {
  #   dir.create(paste("plaatjes/poly_regr/overfit/exp_",experiment_number,sep=""))
  #   setwd(paste("plaatjes/poly_regr/overfit/exp_",experiment_number,sep=""))
  # }
  sd1 = 1
  mean1 = 0
  
  # First, generate sample of order order_sample, but use OLS order order_fit (usually order_fit > order_sample), calculate MSE (will be low)
  if (is.null(slopes)){
    slopes <- runif(order_sample + 1, 0, 1)
  }
  sample <- generate_polynomial_sample(order_p = order_sample, slopes = slopes, size = size, mean = mean1, sd = sd1, x_min = x_min, x_max = x_max, plaatjes = FALSE)
  matrix_x <- sample$matrix_x
  y <- sample$y
  # fit sample using OLS of order order_fit (generaly order_fit >> order_sample in order to demonstrate overfit)
  beta_hat <- ordinair_least_squres_poly(x = matrix_x[,2], y=y, order_p = order_fit, x_min = x_min, x_max = x_max, plaatjes = FALSE)
  
  # plot sample with  regression line
  x <- matrix_x[,2]
  #if (is.null(x_range)) {
  x_range2 <- c(floor(min(x)),ceiling(max(x)))
  #}
  #if (is.null(y_range)) {
  y_range2 <- c(floor(min(y))-1,ceiling(max(y)))
  #}
    
  # write plot to directory
  names <- sprintf("%02d",1:((number_of_test_data_samples+1)*2))
  name <- 1
  #jpeg(file = paste("nozoom_sample_", names[name], ".jpeg", sep=""))
  plot_sample(x = x, y = y, beta_hat = beta_hat, slopes = slopes, x_range = x_range2, y_range = y_range2, sample_line = TRUE)
  #dev.off()
  
  # jpeg(file = paste("outzoom_sample_", names[name], ".jpeg", sep=""))
  if (exists("x_range") && exists("y_range")) {
    print(x_range)
    print("x_range")
    plot_sample(x = x, y = y, beta_hat = beta_hat, slopes = slopes, x_range = x_range, y_range = y_range, sample_line = TRUE)
  }
  # dev.off()
  
  
  # create new matrix_x of order order_fit
  matrix_x = create_matrix_x(matrix_x[,2],order_p = order_fit)
  y_hat <- matrix_x %*% beta_hat
  mse_training_data <- calculate_mse(y = y, y_hat = y_hat)
  # mse_training_data <- mse(matrix_x,y,beta_hat)
  
  # regressie lijn
  y_hat <- seq(0,0,length=dim(matrix_x)[1])
  x_seq <- seq(x_min,x_max,length=dim(matrix_x)[1])
  x_seq_matrix <- create_matrix_x(x_seq, order_p = length(beta_hat))
  for (n in 1:length(y_hat)){
    for (p in 1:(length(beta_hat))){
      #y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq[n]^p
      y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq_matrix[n,p]
    }
  }
  
  # new sample with same slopes, but use same beta hats (test data instead of training data)
  mse_test_data <- c()
  for (k in 1:number_of_test_data_samples){
    sample <- generate_polynomial_sample(order_p = order_sample, slopes = slopes, size = size, mean = mean1, sd = sd1, x_min = x_min, x_max = x_max, plaatjes = FALSE)
    matrix_x <- sample$matrix_x
    y <- sample$y
    x <- matrix_x[,2]
    
    # plot old fit in new data
    # jpeg(file = paste("nozoom_sample_", names[(k+1)], ".jpeg", sep=""))
    # plot_sample(x, y, beta_hat, slopes, x_range = x_range2, y_range2, TRUE)
    print("xrange2")
    print(x_range2)
    plot_sample(x = x, y = y, beta_hat = beta_hat, slopes = slopes, x_range = x_range2, y_range = y_range2, sample_line = TRUE)
    # dev.off()
    if (exists("x_range") && exists("y_range")) {
      # jpeg(file = paste("outzoom_sample_", names[(k+1)], ".jpeg", sep=""))
      plot_sample(x = x, y = y, beta_hat = beta_hat, slopes = slopes, x_range = x_range, y_range = y_range, sample_line = TRUE)
      # dev.off()
    }
    print(paste("xrange",x_range))
    print(paste("yrange",y_range))
    
    # now MSE using beta hats of previous fit, is high
    matrix_x = create_matrix_x(matrix_x[,2],order_p = order_fit)
    y_hat <- matrix_x %*% beta_hat
    mse_test_data[k] <- calculate_mse(y =y, y_hat = y_hat)
    # mse_test_data[k] <- mse(matrix_x,y,beta_hat)
  }
  
  # write data to file
  #write.table(slopes, "plaatjes/poly_regr/overfit/exp1/data.txt",append=TRUE)
  #write.table(paste("beta_hat",beta_hat), "plaatjes/poly_regr/overfit/exp1/data.txt",append=TRUE)
  # write.table(c("slopes"=slopes,"beta_hats"=beta_hat,"mse"=c(mse_training_data, mse_test_data)),"data.txt",append=FALSE)
  # setwd("../../../..")
  return(list("slopes"=slopes,"beta_hats"=beta_hat,"mse"=c(mse_training_data, mse_test_data)))
}

#'
#' Plot data and fit
#' 
plot_sample <- function(x, y, beta_hat, slopes, x_range, y_range, sample_line = FALSE) {
  matrix_x <- create_matrix_x(x,order_p = length(beta_hat))
  plot(x,y, xlim=x_range, ylim=y_range)

  
  # plot regression line based on beta_hat
  plot_length <- 101
  y_hat <- seq(0,0,length=plot_length)
  x_seq <- seq(x_range[1],x_range[2],length=plot_length)
  x_seq_matrix <- create_matrix_x(x_seq, order_p = length(beta_hat))
  
  for (n in 1:length(y_hat)){
    for (p in 1:length(beta_hat)) {
      #y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq[n]^p
      y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq_matrix[n,p]
    }
  }
  #plot(matrix_x[,2],y,ylim=c(0,max(y)+0.1),xlab="x",ylab="y")
  lines(x_seq,y_hat,col="orange")
  
  # plot true sample line
  if (sample_line) {
    y_real <- seq(0,0,length=plot_length)
    for (n in 1:length(y_real)) {
      for (p in 1:length(slopes)) {
        #y_hat[n] <- y_hat[n] + beta_hat[p] * x_seq[n]^p
        y_real[n] <- y_real[n] + slopes[p] * x_seq_matrix[n,p]
      }
    }
    lines(x_seq,y_real,type="l",lty=3)
  }
}

#'
#' Executes cross validation.
#' Data is split into paritions. 
#' Every time one partitons is left out. 
#' The other partitions are moddelled using Bayes with norm prior.
#' Using posterior, predictions for the left out partition are made
#' The is repeated untill every partition is left out once.
#' Result is a vector y_hat with for every y a prediction.
#' y_hat is compared to the observed y by calculating MSE.
#' Return value is MSE and mu_new sigma_new?? To do!
#' 
cross_validation <- function(y, matrix_x, number_of_partitions, alpha){
  # devide sample into partitions
  length_partition <- floor(length(y) / number_of_partitions)
  #print(c("length_partition",length_partition))
  # TO DO: how to handle leftover data points due to rounding?
  #
  partitions_x = list()
  partitions_y = list()
  for (k in 1:number_of_partitions){
    partitions_y[[k]] <- y[(((k-1)*length_partition)+1):(k*length_partition)]
    partitions_x[[k]] <- matrix_x[(((k-1)*length_partition)+1):(k*length_partition),]
  }
  
  # now make training and test data out of partitions
  i <- 1
  y_hats <- seq(0,0,length=length(y))
  
  for (k in 1:number_of_partitions){
    # make sample for training and test data
    
    training_test_x <- make_training_and_test_data_x(number_of_partitions, partitions_x, k)
    training_test_y <- make_training_and_test_data_y(number_of_partitions, partitions_y, k)
    
    training_data_x <- training_test_x$training_data
    test_data_x <- training_test_x$test_data
    training_data_y <- training_test_y$training_data
    test_data_y <- training_test_y$test_data
    
    # now use created training and test data
    posterior <- posterior_normal(matrix_x = training_data_x, y = training_data_y, alpha = alpha, beta_error = 1)
    # To do: what to pick for beta_error? or maybe as variable??
    
    # Take MAP
    beta_hats <- posterior$mu_new
    sigma_new <- posterior$sigma_new
    
    # use beta_hat to estimate y_hats on test data
    
    j <- 1
    stop_i <- i + length_partition -1
    for (i in i:(stop_i)){
      y_hats[i] <- sum(test_data_x[j,] %*% beta_hats)
      j <- j + 1
    }
    i <- i + 1
  }
  #print(c("y_hats",y_hats))
  # To do: delete overtolling zeros at end of y_hats due to rounding
  # calculate mse over all beta hats
  mse <- calculate_mse(y_hats, y)
  
  # calculate also MAP estimate on whole datset
  ## changed on 22-06
  posterior <- posterior_normal(matrix_x = matrix_x, y = y, alpha = alpha, beta_error = 1)
  beta_hats <- posterior$mu_new
  sigma_new <- posterior$sigma_new
  
  return(list("mse"=mse,"mu_new"=beta_hats,"sigma_new"=sigma_new,"y_hats"=y_hats))
}

calculate_mse <- function(y, y_hat){
  mse <- mean((y-y_hat)**2)
  return(mse)
}

#'
#' Help function for cross_valiadation() to split data in test and traing data, based on partitions.
#' 
make_training_and_test_data_x <- function(number_of_partitions, partitions, leave_out){
  
  training_data <- matrix()
  for (k2 in 1:number_of_partitions){
    if (k2 == leave_out) {
      test_data <- partitions[[k2]]
      # test_data_y <- partitions_y[[k2]]
    }
    else {
      # make training data set if now existent
      if (is.na(training_data)[1]){
      training_data <- partitions[[k2]]
      }
      # add data to training data set
      else {
      training_data <- rbind(training_data, partitions[[k2]])
      # training_data_y <- c(training_data_y, partitions_y[[k2]])
      }
    }
    
  }
  return(list("training_data"=training_data,"test_data"=test_data))
}

#'
#' Help function for cross_valiadation() to split data in test and traing data, based on partitions.
#' 
make_training_and_test_data_y <- function(number_of_partitions, partitions, leave_out){
  
  training_data <- c()
  for (k2 in 1:number_of_partitions){
    if (k2 == leave_out) {
      test_data <- partitions[[k2]]
      # test_data_y <- partitions_y[[k2]]
    }
    else {
      # make training data set if now existent
      if (is.null(training_data)){
        training_data <- partitions[[k2]]
      }
      # add data to training data set
      else {
        #training_data <- rbind(training_data, partitions[[k2]])
        training_data <- c(training_data, partitions[[k2]])
      }
    }
    
  }
  return(list("training_data"=training_data,"test_data"=test_data))
}

#'
#'  Gives MAP estimator of normal posterior, using normal prior centred around zero with variance alpha^{-1}
#'  Formulae used from Bishop 3.53 and 3.54
#'  
posterior_normal <- function(matrix_x, y, alpha, beta_error){
  sigma_n_inv <- alpha * diag(ncol(matrix_x)) + beta_error * t(matrix_x) %*% matrix_x
  sigma_n <- solve(sigma_n_inv)
  mu_n <- beta_error * sigma_n %*% t(matrix_x) %*% y
  #mu_new <- sigma_new %*% (sigma_0_inv%*%mu_0 + beta_error*t(matrix_x)%*%y)
  return(list("mu_new"=mu_n,"sigma_new"= sigma_n))
}

#' 
#' Execute the bigg cross validation experiment
#' 
#' 
execute_cross_val <- function(alpha, x_min = -1, x_max = 1, order_sample = 2, order_fit = 2, initial_betas=c(1,1,1), sample_size = c(20), orthog = FALSE, alpha_plot = TRUE, data_plot = TRUE, multiple_samples = FALSE){
  length_alpha <- length(alpha)
  repeat_cross_val <- 1
  number_of_partitions <- 5
  mse_list <- list()
  # mu_list <- list()
  # sigma_list <- list()
  number_of_conditions <- 1
  
  # 2. Fit model using cross validations
  mse_1 <- c()
  sample <- generate_polynomial_sample(order_p = order_sample, slopes = initial_betas, size = sample_size, sd =1, x_min = x_min, x_max = x_max)
  matrix_x_temp <- sample$matrix_x
  y <- sample$y
  x <- matrix_x_temp[,2]
  if (orthog){
    # create orthogmatrix with row of zeros for intercepts (IS THIS CORRECT?? To do!)
    poly1 <- poly(x,order_fit)
    matrix_x_order_fit <- poly1[,]
    matrix_x_order_fit <- cbind(seq(1,1,length=length(x)),matrix_x_order_fit)
    # matrix_x_order_fit <- poly(x,order_fit,simple=TRUE)
    # matrix_x_order_fit <- cbind(seq(0,0,length=n),poly(x,order_fit,simple=TRUE))
  } else {
    matrix_x_order_fit <- create_matrix_x(x = x, order_p = order_fit) 
  }
  for (a in 1:length_alpha){
    mse_2 <- 0
    for (m in 1:repeat_cross_val){
      
      cross_val <- cross_validation(y, matrix_x_order_fit, number_of_partitions, alpha[a])
      mse_2 <- mse_2 + cross_val$mse
      # mu_list[[a]] <- cross_val$mu_new
      # sigma_list[[a]] <- cross_val$sigma_new
      
      # check fit using plot
      check_fit <- FALSE
      y_max <- 1.3 * max(y)
      if (check_fit){
        if (a%%10 == 1 && m == 1){
          # y_max <- 1.3 * max(y)
          plot(matrix_x_order_fit[,2],matrix_x_order_fit%*%cross_val$mu_new,col="red",ylim=c(-3, y_max))
          # legend(0.1,3,paste("mu_list",mu_list[[a]]))
          points(matrix_x_order_fit[,2],y,col="black")
        }
      }
    }
    mse_1[a] <- (mse_2 / repeat_cross_val)
  }
  mse_list[[number_of_conditions]] <- list("betas"=initial_betas,"n"=sample_size,"mses"=mse_1)
  

  if (alpha_plot){

      min_value <- 0.9 * min(mse_list[[number_of_conditions]]$mses)
      max_value <- 1.1 * max(mse_list[[number_of_conditions]]$mses)
      # jpeg(file = paste("plaatjes/cross_val/exp9/","myplot_", names[i], ".jpeg", sep=""))
      plot(alpha,mse_list[[number_of_conditions]]$mses,ylim=c(min_value,max_value),log="",xlab=expression(alpha),ylab="MSE")
      # lines(seq(0,max(alpha),length=repeat_cross_val),rep(round(mean(mse_list[[i]]$mses),3),length=repeat_cross_val),lty=3,lwd=3,col="red")
      
      # legend(0.1*max(alpha),max_value-0.05,paste(c("betas",mse_list[[i]]$betas)))
      # legend(0.3*max(alpha),max_value-0.05,paste("n",mse_list[[i]]$n))
      # legend(0.5*max(alpha),max_value-0.05,paste(c("mean",round(mean(mse_list[[i]]$mses),3))))
      # legend(0.7*max(alpha),max_value-0.05,"fit 2, true 2, folds 5")
      # dev.off()

  }
  #if (data_plot && !multiple_samples){
  if (!multiple_samples){
    # var for plot
    x_seq <- seq(x_min,x_max,length=length(x))
    if (orthog){
      # no intercepts
      x_poly <- poly(x,order_fit)
      
      # print(poly(x_seq, order_fit, coefs=attr(x_poly,"coefs"))[,])
      matrix_x_seq <- x_poly[,]
      # print(poly(x_seq, order_p, coefs=attr(x_poly,"coefs"))[,])
      # matrix_x_seq <- poly(x_seq, order_p, coefs=attr(x_poly,"coefs"))[,]
      # matrix_x_seq <-create_matrix_x(x = x_seq, order_p = order_fit)#[,2:ncol(matrix_x_seq)]
      # matrix_x_seq <- matrix_x_seq[,2:ncol(matrix_x_seq)]
    } else {
      matrix_x_seq <-create_matrix_x(x = x_seq, order_p = order_fit)#[,2:3]
    }
    
    
    # find alpha*
    min_index <- which(mse_list[[number_of_conditions]]$mses==min(mse_list[[number_of_conditions]]$mses))
    alpha_star <- alpha[min_index]

    cross_val_alpha_star <- cross_validation(y, matrix_x_order_fit, number_of_partitions, alpha_star)
    #  compare to low alha
    cross_val_alpha_laag <- cross_validation(y, matrix_x_order_fit, number_of_partitions, alpha[1])
    # compare to high alpha
    cross_val_alpha_hoog <- cross_validation(y, matrix_x_order_fit, number_of_partitions, alpha[length(alpha)])
    
    if (data_plot){
      
      # alpha_star_post <- posterior_normal(matrix_x_order_fit, y, alpha =alpha_star, beta_error =1)
      if (orthog) {
        print(matrix_x_seq)
        print(seq(1,1,length=length(x)))
        matrix_x_seq_temp <- matrix_x_seq[order(matrix_x_seq[,1]),]
        matrix_x_seq <- cbind(seq(1,1,length=length(x)),matrix_x_seq_temp)
        print(matrix_x_seq)
        
      }
       
      plot(matrix_x_order_fit[,2],y,col="black",ylim=c(-8,y_max),xlab = "x",ylab="y")
      alpha_star_post <- posterior_normal(matrix_x_order_fit, y, alpha =alpha_star, beta_error =1)
      lines(matrix_x_seq[,2],matrix_x_seq%*%alpha_star_post$mu_new,col="green")
      points(matrix_x_order_fit[,2],matrix_x_order_fit%*%alpha_star_post$mu_new,col="green")
      
      print("alpha star ")
      print(alpha_star)
      print(alpha_star_post)
      
      # alpha low
      alpha_laag_post <- posterior_normal(matrix_x_order_fit, y, alpha =alpha[1], beta_error =1)
      lines(matrix_x_seq[,2],matrix_x_seq%*%alpha_laag_post$mu_new,col="red", type = "l")
      points(matrix_x_order_fit[,2],matrix_x_order_fit%*%alpha_laag_post$mu_new,col="red")
      
      
      # points(matrix_x_order_fit[,2],matrix_x_order_fit%*%cross_val_alpha_hoog$mu_new,col="orange")
      #lines(x_seq,matrix_x_seq%*%cross_val_alpha_hoog$mu_new,col="orange")
      
      alpha_hoog_post <- posterior_normal(matrix_x_order_fit, y, alpha =alpha[length(alpha)], beta_error =1)
      lines(matrix_x_seq[,2],matrix_x_seq%*%alpha_hoog_post$mu_new,col="orange")
      points(matrix_x_order_fit[,2],matrix_x_order_fit%*%alpha_hoog_post$mu_new,col="orange")
      
      # # quick lm
      # dataframe_alphahoog <- data.frame("y_hats"=cross_val_alpha_hoog$y_hats, "x"=x, "x**2"=x**2)
      # model_alphahoog <- lm(y ~ x + x..2,dataframe_alphahoog)
      # print(model_alphahoog)
      # print(dataframe_alphahoog)
      # lines(x_seq,model_alphahoog$coef[[1]] + model_alphahoog$coef[[2]] * x_seq + model_alphahoog$coef[[3]] * x_seq**2,col="orange")
      
      # legend
      legend(x="topleft",legend=c("Ware data","Voorspelling alpha*","Voorspelling lage alpha","Voorspelling hoge alpha" ), col=c("black","green","red","orange"),lwd=1,pch=1,cex=0.75)
    }

  }
  if (!multiple_samples){
    mse_list[[number_of_conditions]][["alpha_star"]] <- alpha_star
    mse_list[[number_of_conditions]][["alpha_star_cross_val"]] <- cross_val_alpha_star
    mse_list[[number_of_conditions]][["alpha_laag"]] <- alpha[1]
    mse_list[[number_of_conditions]][["alpha_laag_cross_val"]] <- cross_val_alpha_laag
    mse_list[[number_of_conditions]][["alpha_hoog"]] <- alpha[length(alpha)]
    mse_list[[number_of_conditions]][["alpha_hoog_cross_val"]] <- cross_val_alpha_hoog
  }

  return(mse_list)
}

#'
#' Calculates out of back MSE
#' Using beta_hats a new random sample is generated.
#' This random sample is compared to the true prediction (without the error term)
#' 
out_of_back_mse <- function(beta_hat, beta_true, x_min, x_max, n = 100){
  
  # variables
  order_p <- length(beta_true) - 1
  order_fit <- length(beta_hat) - 1
  
  # generate sample
  sample1 <- generate_polynomial_sample(order_p = order_p, slopes=beta_true, size=n, mean=0, sd=1, x_min=x_min,x_max=x_max, plaatjes=FALSE)
  matrix_x1 <- sample1$matrix_x
  y1 <- sample1$y
  
  matrix_x1_order_fit <- create_matrix_x(x = matrix_x1[,2], order_p = order_fit)
  
  # calculate y_hat
  # x_seq <- seq(x_min,x_max,length=n)
  # x_seq_matrix <- create_matrix_x(x = x_seq, order_p = order_fit)
  # y_hat <- x_seq_matrix %*% beta_hat
  y_hat <- matrix_x1_order_fit %*% beta_hat
  
  # calculate mse
  mse <- calculate_mse(y = y1, y_hat = y_hat)
  
  return(mse)
}
