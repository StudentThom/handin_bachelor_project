###############################################################################
#' helpers_CV_bayesglm.R 													  #
#' helpers file for basic and polynomial bayesian regression using bayesgl	  #
#' i.e. T-priors															  #
#' Created by StudentThom, spring 2018										  #
###############################################################################
#'
#' In this file:
#'
#' execute_cross_val_tprior(df, x_min = -1, x_max = 1, order_sample = 2, order_fit = 2,
#'
#' cross_validation_tprior(y, matrix_x, number_of_partitions, df)
#'

library("arm")

#
execute_cross_val_tprior <- function(df, x_min = -1, x_max = 1, order_sample = 2, order_fit = 2, initial_betas=(c(1,1,1)), sample_size = c(20), orthog = FALSE, df_plot = TRUE, data_plot = TRUE, multiple_samples = FALSE){
  length_df <- length(df)
  repeat_cross_val <- 1
  number_of_partitions <- 5
  number_of_conditions <- 1
  mse_list <- list()


  # 2. Fit model using cross validations
  mse_1 <- c()
  sample <- generate_polynomial_sample(order_p = order_sample, slopes = initial_betas, size = sample_size, sd =1, x_min = x_min, x_max = x_max)
  matrix_x_temp <- sample$matrix_x
  y <- sample$y
  x <- matrix_x_temp[,2]
  y_max <- max(y) * 1.1
  if (orthog){
    # create orthogmatrix with row of zeros for intercepts (IS THIS CORRECT?? To do!)
    poly1 <- poly(x,order_fit)
    matrix_x_order_fit <- poly1[,]
    # matrix_x_order_fit <- poly(x,order_fit,simple=TRUE)
    matrix_x_order_fit <- cbind(seq(0,0,length=n),poly(x,order_fit,simple=TRUE))
  } else {
    matrix_x_order_fit <- create_matrix_x(x = x, order_p = order_fit) 
  }
  for (a in 1:length_df){
    mse_2 <- 0

    cross_val <- cross_validation_tprior(y, matrix_x_order_fit, number_of_partitions, df[a],order_fit = order_fit)
    mse_2 <- mse_2 + cross_val$mse
    # mu_list[[a]] <- cross_val$mu_new
    # sigma_list[[a]] <- cross_val$sigma_new

    mse_1[a] <- (mse_2 / repeat_cross_val)
  }
  mse_list[[number_of_conditions]] <- list("betas"=initial_betas,"n"=sample_size,"mses"=mse_1)
  
  
  if (df_plot){
    
    min_value <- 0.9 * min(mse_list[[number_of_conditions]]$mses)
    max_value <- 1.1 * max(mse_list[[number_of_conditions]]$mses)
    # jpeg(file = paste("plaatjes/cross_val/exp9/","myplot_", names[i], ".jpeg", sep=""))
    plot(df,mse_list[[number_of_conditions]]$mses,ylim=c(min_value,max_value),log="",xlab=expression(df),ylab="MSE")
    # lines(seq(0,max(df),length=repeat_cross_val),rep(round(mean(mse_list[[i]]$mses),3),length=repeat_cross_val),lty=3,lwd=3,col="red")
    # regression line for when multiple samples where used
    if (multiple_samples){
      mse_regr <- lm(mse_list[[1]]$mses ~ df)
      mse_intercept <- mse_regr$coef[[1]]
      mse_slope <- mse_regr$coef[[2]]
      lines(df,mse_intercept + mse_slope * df,col="red")
    }
    # legend(0.1*max(alpha),max_value-0.05,paste(c("betas",mse_list[[i]]$betas)))
    # legend(0.3*max(alpha),max_value-0.05,paste("n",mse_list[[i]]$n))
    # legend(0.5*max(alpha),max_value-0.05,paste(c("mean",round(mean(mse_list[[i]]$mses),3))))
    # legend(0.7*max(alpha),max_value-0.05,"fit 2, true 2, folds 5")
    # dev.off()
    
  }
  #if (data_plot && !multiple_samples){
  if (!multiple_samples){
    # var for plot
    x_seq <- seq(x_min,x_max,length=101)
    if (orthog){
      # no intercepts
      matrix_x_seq <-create_matrix_x(x = x_seq, order_p = order_fit)#[,2:ncol(matrix_x_seq)]
      matrix_x_seq <- matrix_x_seq[,2:ncol(matrix_x_seq)]
    } else {
      matrix_x_seq <-create_matrix_x(x = x_seq, order_p = order_fit)#[,2:3]
    }
    
    
    # find df*
    min_index <- which(mse_list[[number_of_conditions]]$mses==min(mse_list[[number_of_conditions]]$mses))
    df_star <- df[min_index]
    
    cross_val_df_star <- cross_validation_tprior(y, matrix_x_order_fit, number_of_partitions, df_star, order_fit = order_fit)
    #  compare to low alha
    cross_val_df_laag <- cross_validation_tprior(y, matrix_x_order_fit, number_of_partitions, df[1], order_fit = order_fit)
    # compare to high df
    cross_val_df_hoog <- cross_validation_tprior(y, matrix_x_order_fit, number_of_partitions, df[length(df)], order_fit = order_fit)
    
    if (data_plot){
      plot(matrix_x_order_fit[,2],y,col="black",ylim=c(-3,1.2*y_max),xlab = "x",ylab="y")
      # alpha *
      # points(matrix_x_order_fit[,2],matrix_x_order_fit%*%cross_val_alpha_star$mu_new,col="green")
      lines(x_seq,matrix_x_seq%*%cross_val_df_star$mu_new,col="green")
      # legend(0.1,3,paste("mu_list",mu_list[[a]]))
      lines(x_seq,matrix_x_seq%*%cross_val_df_laag$mu_new,col="red")
      
      
      # points(matrix_x_order_fit[,2],matrix_x_order_fit%*%cross_val_df_hoog$mu_new,col="orange")
      lines(x_seq,matrix_x_seq%*%cross_val_df_hoog$mu_new,col="orange")
      # legend
      legend(x="topleft",legend=c("Ware data","Voorspelling df*","Voorspelling lage df","Voorspelling hoge df" ), col=c("black","green","red","orange"),lwd=1,pch=1,cex=0.75)
    }
    
  }
  if (!multiple_samples){
    mse_list[[number_of_conditions]][["df_star"]] <- df_star
    mse_list[[number_of_conditions]][["df_star_cross_val"]] <- cross_val_df_star
    mse_list[[number_of_conditions]][["df_laag"]] <- df[1]
    mse_list[[number_of_conditions]][["df_laag_cross_val"]] <- cross_val_df_laag
    mse_list[[number_of_conditions]][["df_hoog"]] <- df[length(df)]
    mse_list[[number_of_conditions]][["df_hoog_cross_val"]] <- cross_val_df_hoog
  }

  return(mse_list)
}


#################################
#
cross_validation_tprior <- function(y, matrix_x, number_of_partitions, df, order_fit){
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
    training_data_dataframe <- data.frame(y=training_data_y, training_data_x)
    # # print(training_data_dataframe)
    # posterior <- bayesglm(formula = y ~ X2 + X3 + X4 + X5 +X6, data = training_data_dataframe,prior.df=df)$coef
    # # print(posterior)
    # beta_hats <- c(posterior[[1]],posterior[[2]],posterior[[3]],posterior[[4]],posterior[[5]], posterior[[6]])
    
    if (order_fit == 2){
      posterior <- bayesglm(formula = y ~ X2 + X3, data = training_data_dataframe,prior.df=df)$coef
      beta_hats <- c(posterior[[1]],posterior[[2]],posterior[[3]])
    }
    else if (order_fit == 5){
      posterior <- bayesglm(formula = y ~ X2 + X3 + X4 + X5 +X6, data = training_data_dataframe,prior.df=df)$coef
      beta_hats <- c(posterior[[1]],posterior[[2]],posterior[[3]],posterior[[4]],posterior[[5]], posterior[[6]])
    }
    
    # use beta_hat to estimate y_hats on test data
    j <- 1
    stop_i <- i + length_partition -1
    for (i in i:(stop_i)){
      # print("test data")
      # print(test_data_x)
      #print(beta_hats)
      y_hats[i] <- sum(test_data_x[j,] %*% beta_hats)
      j <- j + 1
    }
    i <- i + 1
  }
  #print(c("y_hats",y_hats))
  # To do: delete overtolling zeros at end of y_hats due to rounding
  # calculate mse over all beta hats
  mse <- calculate_mse(y_hats, y)
  
  # to do
  if (order_fit == 2){
    posterior <- bayesglm(formula = y ~ X2 + X3, data = training_data_dataframe,prior.df=df)$coef
    beta_hats <- c(posterior[[1]],posterior[[2]],posterior[[3]])
  }
  else if (order_fit == 5){
    posterior <- bayesglm(formula = y ~ X2 + X3 + X4 + X5 +X6, data = training_data_dataframe,prior.df=df)$coef
    beta_hats <- c(posterior[[1]],posterior[[2]],posterior[[3]],posterior[[4]],posterior[[5]], posterior[[6]])
  }

  
  return(list("mse"=mse,"mu_new"=beta_hats))
}

