# File Containing all functions that deal with Nomograms
# By: Mohammed Janahi
#
#

LEVELS <- c(2.5 , 5, 10 , 25 , 50 , 75 , 90 , 95 , 97.5)

# Main Nomogram Function.
# Function to perform the sliding window analysis
#   After thorough investigation of the Nomogram Paper, I can say with some degree of confidence that the way Nomograms
#   are generated is: 
#   Given a table of samples, with Ages and HV columns,  ( lets assume they are 4000 samples )
#   At every 1% of samples start creating bins/windows. (0% to 100% so bin start indices will go from 1 to 3960 )
#   For every bin/window, select the set bin size number of samples from the start index.
#   Bin size is given as a percent of the total number of samples. (so, if the percent is 10, then samples per bin would be 400)
# INPUTS:
#   required:
#       - ukb: Table of samples with at age, left hippocampal volume, and right hippocampal volume 
#   optional:
#       - percent_samples_per_bin: percent of samples to be included in each bin, as a percent of total number of samples (0-100)
#       - kernel_width : width of the kernel to be used for Gaussian smoothing
#       - age_column : name of the column with the ages
#       - hv_left_column : name of the column with the left hippocampal volume
#       - hv_right_column : name of the column with the right hippocampal volume
# OUTPUT:
#   table of results for each bin/window. one bin/window per row and columns are:
#    min_age : minimum age of samples included in this bin/window
#    max_age : maximum age of samples included in this bin/window
#    mid_age : median age of samples included in this bin/window

#    mean_left : mean of left hippocampal volume of samples included in this bin/window
#    std_left : standard deviation of left hippocampal volume of samples included in this bin/window
#    min_left :  minimum left hippocampal volume of samples included in this bin/window
#    q2.5_left : 
#    q5_left : 
#    ...
#    q97.5_left :  
#    max_left : maximum left hippocampal volume of samples included in this bin/window

#    mean_right : same as before, but for right HV
#    ...
#    max_right : same as before, but for right HV

#    mean : same as before, but for mean HV
#    ...
#    max : same as before, but for mean HV

# ASSUMPTIONS:
#
WINDOW_ANALYSIS <- function( ukb , percent_samples_per_bin=10 , kernel_width = 20,
                             age_column="AGE_Latest" ,
                             hv_mean_column="clean_hv_bilateral" ,
                             hv_left_column="clean_hv_left", 
                             hv_right_column="clean_hv_right"){
  
  # sort samples by age
  ukb <- ukb[ order(ukb[,age_column]), ]
  
  NUM_SAMPLES <- nrow(ukb)
  
  # calculate the size of each bin
  samples_per_bin <- NUM_SAMPLES / percent_samples_per_bin
  
  # 1% of samples will determine the start  index of each bin
  one_percent_of_samples <- NUM_SAMPLES/100
  
  # initialize empty output table
  bins <- INITIALIZE_BINS()
  
  # loop from beginning, increment by 1% of samples each time, stop before the size of the last bin exceeds num of samples
  for (i in seq( 1 , (NUM_SAMPLES - samples_per_bin) , one_percent_of_samples) ){
    # get samples on this window
    window_of_samples <- ukb[i:(i+samples_per_bin) , c(age_column, hv_left_column , hv_right_column , hv_mean_column)]
    
    win_ages      <- window_of_samples[ , c(age_column) ]
    win_min_age   <- min(win_ages  , na.rm = TRUE)
    win_max_age   <- max(win_ages  , na.rm = TRUE)
    win_mean_age  <- mean(win_ages , na.rm = TRUE)
    
    win_hv_left   <- window_of_samples[,hv_left_column]
    win_hv_right  <- window_of_samples[,hv_right_column]
    win_hv  <- window_of_samples[,hv_mean_column]
    
    win_mean_hv_left   <- mean(win_hv_left, na.rm = TRUE) 
    win_std_hv_left    <- sd(win_hv_left , na.rm = TRUE)
    win_min_hv_left    <- min(win_hv_left , na.rm = TRUE)
    win_max_hv_left    <- max(win_hv_left , na.rm = TRUE)
    win_quantiles_left <- quantile( win_hv_left , probs = LEVELS/100 , na.rm = TRUE) 
    
    win_mean_hv_right   <- mean(win_hv_right, na.rm = TRUE) 
    win_std_hv_right    <- sd(win_hv_right , na.rm = TRUE)
    win_min_hv_right    <- min(win_hv_right , na.rm = TRUE)
    win_max_hv_right    <- max(win_hv_right , na.rm = TRUE)
    win_quantiles_right <- quantile( win_hv_right , probs = LEVELS/100 , na.rm = TRUE) 

    win_mean_hv   <- mean(win_hv, na.rm = TRUE) 
    win_std_hv    <- sd(win_hv , na.rm = TRUE)
    win_min_hv    <- min(win_hv , na.rm = TRUE)
    win_max_hv    <- max(win_hv , na.rm = TRUE)
    win_quantiles <- quantile( win_hv , probs = LEVELS/100 , na.rm = TRUE) 
    
    row <- c( win_min_age , win_max_age , win_mean_age ,
              win_mean_hv_left , win_std_hv_left ,
              win_min_hv_left , win_quantiles_left , win_max_hv_left ,
              win_mean_hv_right , win_std_hv_right,
              win_min_hv_right , win_quantiles_right , win_max_hv_right,
              win_mean_hv, win_std_hv,
              win_min_hv , win_quantiles , win_max_hv )

    bins[nrow(bins)+1,] <- row
  }
  
  # Gaussian smoothing
  for (col in names(bins)){
    bins[ , col ] <- smth( bins[, col ],  method = 'gaussian', window=kernel_width)
  }
  
  # return the non-na rows
  return(bins[!is.na(bins$mid_age),])
}

#
#
#
# IF YOU GIVE IT LEFT AND RIGHT COLUMNS, IT WILL RETRAIN FOR THOSE AS WELL.
GPR_ANALYSIS <- function(ukb , XX=NA , age_column="AGE_Latest" ,
                         hv_column="clean_hv_bilateral" , hv_left_column=NA, hv_right_column=NA , kernel_width=20){

  if(is.na(XX[1])){
    XX <- matrix(seq( 30, 100, length = 70*4) , ncol=1)
  }
  
  # initialize the output table
  bins <- INITIALIZE_BINS(nrow(XX))
  
  bins$mid_age <- c(XX)
  
  for( hem in na.omit(c(hv_column , hv_left_column , hv_right_column )) ){
    
    suf <- ifelse(identical(hem, hv_left_column) , "_left" , (ifelse( identical(hem, hv_right_column) , "_right" , "" )) )
    
    y = ukb[ !is.na(ukb[,hem]) , hem]
    X = ukb[ !is.na(ukb[,hem]) , age_column]
    
    # experimental bit adding regression line as samples 
    #a <- lm(y~X)
    #alpha <- a$coefficients[1]
    #beta <- a$coefficients[2]
    #X_lm <- seq( 30, 100, length = 70*4)
    #y_lm <- alpha + (X_lm * beta)
    #X <- c(X , X_lm)
    #y <- c(y , y_lm)
    
    g2 <- garg(list(mle = TRUE , start=mean(y)), y )
    d2 <- darg(list(mle = TRUE, start=mean(distance(X)) ), matrix(X , ncol=1) )
    
    gp  <- myGPsep(matrix(X , ncol=1) ,y , d=d2$start , g= g2$start , dK=TRUE) 
    mle <- mleGPsep(  gp$gpsepi , param="both" , tmin=c(d2$min,g2$min) , tmax=c(d2$max, g2$max) )
    p   <- predGPsep( gp$gpsepi , XX )
    
    deleteGPsep(gp$gpsepi)
    
    bins[ , paste("mean",suf,sep="") ] <- p$mean
    bins[ , paste("std",suf,sep="") ] <- sqrt(diag(p$Sigma))
    for(lvl in LEVELS)
      bins[ , paste( "q",lvl,suf,sep="")]  <- qnorm(lvl/100, mean = p$mean, sd = sqrt(diag(p$Sigma)))
    
    bins[ , paste("min",suf,sep="") ] <- -Inf
    bins[ , paste("max",suf,sep="") ] <- Inf
  }
  
  # Gaussian smoothing
  for (col in names(bins)){
    if(!is.na(bins[ 1, col ]))
      bins[ , col ] <- smth( bins[, col ],  method = 'gaussian', window=kernel_width)
  }
  
  # return the non-na rows
  return(bins[!is.na(bins$mid_age),])
}


# my edit of the newGPsep function from the laGP package.
# edited to return the full GP and not just the index.
myGPsep <- function (X, Z, d, g, dK = FALSE) {
  n <- nrow(X)
  m <- ncol(X)
  if (is.null(n)) 
    stop("X must be a matrix")
  if (length(Z) != n) 
    stop("must have nrow(X) = length(Z)")
  if (length(d) == 1) 
    d <- rep(d, m)
  else if (length(d) != m) 
    stop("must have length(d) = ncol(X)")
  out <- .C("newGPsep_R", m = as.integer(m), n = as.integer(n), 
            X = as.double(t(X)), Z = as.double(Z), d = as.double(d), 
            g = as.double(g), dK = as.integer(dK), gpsepi = integer(1), 
            PACKAGE = "laGP")
  return(out)
}


#
#
#
GPR_MODEL <- function(ukb , x_cols=c("AGE_Latest") , y_col="clean_hv_bilateral" ){

  y = ukb[ !is.na(ukb[,y_col]) , y_col]
  X = ukb[ !is.na(ukb[,y_col]) , x_cols]

  g2 <- garg(list(mle = TRUE), y)
  d2 <- darg(list(mle = TRUE, max = 100), as.matrix(X) )
  
  gp <- myGPsep(as.matrix(X) ,y , d=d2$start , g= g2$start , dK=TRUE) 
  mle <- mleGPsep(gp$gpsepi , param="both" , tmin=c(d2$min,g2$min) , tmax=c(d2$max, g2$max) )
  
  return(gp)
}

GPR_MODEL_TO_BINS <- function(gp , XX){
  
  p <- predGPsep( gp$gpsepi , XX )
  
  if(is.na(XX[1])){
    XX <- matrix(seq( 30, 100, length = 70*4) , ncol=1)
  }
  
  # initialize the output table
  bins <- INITIALIZE_BINS(nrow(XX))
  
  bins$mid_age <- c(XX)
  bins[ , paste("mean",suf,sep="") ] <- p$mean
  bins[ , paste("std",suf,sep="") ] <- sqrt(diag(p$Sigma))
  for(lvl in LEVELS)
    bins[ , paste( "q",lvl,suf,sep="")]  <- qnorm(lvl/100, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  
  bins[ , paste("min",suf,sep="") ] <- -Inf
  bins[ , paste("max",suf,sep="") ] <- Inf
}

GPR_ANALYSIS_TEST_MF <- function(ukb , XX=NA , age_column="AGE_Latest" , 
                              hv_left_column="clean_hv_left", hv_right_column="clean_hv_right"){
  if(is.na(XX[1])){
    L = 70*4
    XX_male <- c( seq( 30, 100, length = L) ,  (1:L*0 + 1) )
    dim(XX_male) <- c(L,2)
    XX_female <- c( seq( 30, 100, length = L) ,  (1:L*0 - 1) )
    dim(XX_female) <- c(L,2)
  }
  
  # initialize the output table
  bins_male <- INITIALIZE_BINS(nrow(XX))
  
  # initialize the output table
  bins_female <- INITIALIZE_BINS(nrow(XX))
  
  
  X = ukb[ , c(age_column,"PRS_TH_1_Levels")]
  y = ukb[ , hv_left_column]
  
  X = X[!is.na(y) , ]
  y = y[!is.na(y)]
  
  bins_male$mid_age <- XX_male[,1]
  bins_female$mid_age <- XX_male[,1]
  
  gpi <- newGPsep(as.matrix(X) , y , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  
  p <- predGPsep( gpi , as.matrix(XX_male) )
  
  bins_male$mean_left <- p$mean
  bins_male$std_left <- sqrt(diag(p$Sigma))
  bins_male$l.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$l.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$min_left <- -Inf
  bins_male$max_left <- Inf
  
  p <- predGPsep( gpi , as.matrix(XX_female) )
  
  bins_female$mean_left <- p$mean
  bins_female$std_left <- sqrt(diag(p$Sigma))
  bins_female$l.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$l.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$min_left <- -Inf
  bins_female$max_left <- Inf
  
  deleteGPsep(gpi)
  
  X = ukb[ , c(age_column,"PRS_TH_1_Levels")]
  y = ukb[ , hv_right_column]
  
  X = X[!is.na(y) , ]
  y = y[!is.na(y)]
  
  
  gpi <- newGPsep(as.matrix(X) ,y , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  
  p <- predGPsep( gpi , as.matrix(XX_male) )
  
  bins_male$mean_right <- p$mean
  bins_male$std_right <- sqrt(diag(p$Sigma))
  bins_male$r.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$r.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_male$min_right <- -Inf
  bins_male$max_right <- Inf
  
  p <- predGPsep( gpi , as.matrix(XX_female) )
  
  bins_female$mean_right <- p$mean
  bins_female$std_right <- sqrt(diag(p$Sigma))
  bins_female$r.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$r.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins_female$min_right <- -Inf
  bins_female$max_right <- Inf
  
  deleteGPsep(gpi)
  
  return( list(bins_high = bins_male,bins_low = bins_female) )
}

GPR_ANALYSIS_TEST_LR <- function(ukb , XX=NA , age_column="AGE_Latest" , hv_column="clean_hv"){

    L = nrow(ukb)/2
  if(is.na(XX[1])){
    XX_R <- c( seq( 30, 100, length = L) ,  (seq( 30, 100, length = L)*0 + 1) )
    dim(XX_R) <- c(L,2)
      XX_L <- c( seq( 30, 100, length = L) ,  (seq( 30, 100, length = L)*0 ) )
    dim(XX_L) <- c(L,2)
  }
  
  # initialize the output table
  bins <- INITIALIZE_BINS( L )

  X = ukb[ , c(age_column,"L/R")]
  y = ukb[ , hv_column]
  
  X = X[!is.na(y) , ]
  y = y[!is.na(y)]
  
  bins$mid_age <- XX_L[,1]

  gpi <- newGPsep(as.matrix(X) , y , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  
  p <- predGPsep( gpi , as.matrix(XX_L) )
  
  bins$mean_left <- p$mean
  bins$std_left <- sqrt(diag(p$Sigma))
  bins$l.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$l.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$min_left <- -Inf
  bins$max_left <- Inf
  
  p <- predGPsep( gpi , as.matrix(XX_R) )
  
  bins$mean_right <- p$mean
  bins$std_right <- sqrt(diag(p$Sigma))
  bins$r.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$r.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  bins$min_right <- -Inf
  bins$max_right <- Inf
  
  deleteGPsep(gpi)
  
  return(bins)
}


GPR_ANALYSIS_2D <- function(){
  
  library(laGP)
  X = ukb_img_ex_outliers_male[ , c("AGE_Latest","PRS_TH_1")]
  y = ukb_img_ex_outliers_male$clean_hv_left
  
  X = X[!is.na(y) , ]
  y = y[!is.na(y)]
  
  X = X[!is.na(X$PRS_TH_1) , ]
  y = y[!is.na(X$PRS_TH_1)]
  
  XX <- as.matrix( expand.grid( seq( 30, 90, length = 100) , seq( -10, 10, length = 100) ))
  
  gpi <- newGPsep( as.matrix(X[1:1000 , ] ) , y[1:1000] , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  p <- predGPsep( gpi , XX )
  deleteGPsep(gpi)
  
  YY <- rmvnorm( 100 , p$mean , p$Sigma)
  
  matplot(XX , t(YY) , type = "l" , col="grey" , lty = 1)
  lines(XX, p$mean, col = 2, lwd = 2)
  points(X,y , col=rgb(0,1,0,0.1) )
  q1 <- qnorm(0.05, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  q2 <- qnorm(0.95, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  lines(XX, q1, col = 4, lty = 2, lwd = 2 )
  lines(XX, q2, col = 4, lty = 2, lwd = 2 )
  
  
}


GPR_WITH_PLOT <- function(){
  X = ukb_img_ex_outliers_male$AGE_Latest
  y = ukb_img_ex_outliers_male$clean_hv_left
  
  X = matrix(X[!is.na(y)] , ncol=1)
  y = y[!is.na(y)]
  
  
  XX <- matrix(seq( 30, 90, length = 100),ncol=1)
  
  gpi <- newGPsep(matrix(X[1:100] , ncol=1) ,y[1:100] , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  p <- predGPsep( gpi , XX )
  deleteGPsep(gpi)
  
  YY <- rmvnorm( 100 , p$mean , p$Sigma)
  
  matplot(XX , t(YY) , type = "l" , col="grey" , lty = 1)
  lines(XX, p$mean, col = 2, lwd = 2)
  points(X,y , col=rgb(0,1,0,0.1) )
  q1 <- qnorm(0.05, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  q2 <- qnorm(0.95, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  lines(XX, q1, col = 4, lty = 2, lwd = 2 )
  lines(XX, q2, col = 4, lty = 2, lwd = 2 )
  
}


INITIALIZE_BINS <- function( length=0 ){

  bins <- data.frame(matrix(ncol = 42, nrow = length))
  names(bins) <- c("min_age" , "max_age" , "mid_age" , "mean_left" , "std_left" , "min_left" ,
                      paste("q" , LEVELS , "_left"  ,sep=""), "max_left","mean_right","std_right","min_right", 
                      paste("q" , LEVELS , "_right" ,sep=""), "max_right","mean","std","min", 
                      paste("q" , LEVELS , sep=""), "max" )
  return(bins)
}

#Assumptions:
# - Bins have the exact same mid_age column entries (or very close)
# - Bins are constructed as defined in the ANALYSIS functions (they have columns l.q* and r.q* )
NOMOGRAM_DIFF <- function( bins_1 , bins_2 , hem="bilateral" ){
  
  suf <- ifelse( hem=="bilateral" , "", paste("_",hem,sep=""))
  cols <- paste( "q",LEVELS,suf ,sep="")

  mean_error <- mean(colMeans(abs( bins_1[ , cols] - bins_2[ , cols])))
  
  return(mean_error)
}

NOMOGRAM_DIFF_INTERPOLATE <- function( bins_1 , bins_2 , hem="bilateral" , age_range = NA ){

  error = 0
  
  suf <- ifelse( hem=="bilateral" , "", paste("_",hem,sep=""))
  cols <- paste( "q",LEVELS,suf ,sep="")
  
  if( is.na(age_range[1]) ){
    age_range_min <- max(min(bins_1$mid_age),min(bins_2$mid_age))
    age_range_max <- min(max(bins_1$mid_age),max(bins_2$mid_age))
  }
  else{
    age_range_min = age_range[1]
    age_range_max = age_range[2]
  }
  for ( col_name in cols ){
    
    for( i in 1:nrow(bins_1) ){
      
      age = bins_1[i , "mid_age"]
      volume_1 = bins_1[ i , col_name ]
      if(  age >= age_range_min  & age <= age_range_max ){
        min_ind <- which.min( abs(age - bins_2$mid_age) )
        
        if( bins_2[ min_ind , "mid_age"] <= age ){
          under_age <- bins_2[ min_ind , "mid_age"]
          over_age <- bins_2[ min_ind+1 , "mid_age"]
          
          under_volume <- bins_2[ min_ind , col_name]
          over_volume <- bins_2[ min_ind+1 , col_name]
        }
        else{
          over_age <- bins_2[ min_ind , "mid_age"]
          under_age <- bins_2[ min_ind-1 , "mid_age"]
          
          over_volume <- bins_2[ min_ind , col_name]
          under_volume <- bins_2[ min_ind-1 , col_name]
        }
        
        volume_2 <-  (( (age - under_age) / (over_age - under_age) ) * (over_volume - under_volume) ) + under_volume
        if(!is.na(abs( volume_1 - volume_2)))
          error <- error + abs( volume_1 - volume_2)
      }
    }
    
  }
  
  error <- error / (length(cols) * nrow(bins_1) )
    #error = error + sum( abs( bins_1[ , col_name] - bins_2[ , col_name] ) , na.rm = TRUE )
  max_age_range <- max(max(bins_1$mid_age),max(bins_2$mid_age)) - min(min(bins_1$mid_age),min(bins_2$mid_age))
  min_age_range <- (age_range_max - age_range_min)
  overlap_percent <- min_age_range / max_age_range 
  #print(paste("AGE RANGE OVERLAP:" , overlap_percent))
  return(error)
}

# Function to plot nomograms given a bins object
#
PLOT_NOMOGRAM <- function( bins , hem = "bilateral" , title = " " , xlim = NA , ylim = NA ){
  
  centers <- bins$mid_age
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )

  mn <- bins[ , paste("mean",suf,sep="") ]
  q2.5 <- bins[ , paste("q2.5",suf,sep="")]
  q5 <- bins[ , paste("q5",suf,sep="")]
  q10 <- bins[ , paste("q10",suf,sep="")]
  q25 <- bins[ , paste("q25",suf,sep="")]
  q50 <- bins[ , paste("q50",suf,sep="")]
  q75 <- bins[ , paste("q75",suf,sep="")]
  q90 <- bins[ , paste("q90",suf,sep="")]
  q95 <- bins[ , paste("q95",suf,sep="")]
  q97.5 <- bins[ , paste("q97.5",suf,sep="")]
  

  if( is.na(xlim[1])){
    min_x <- floor(head(centers , 1 )) - 2
    max_x <- floor(tail(centers , 1 )) + 2
  }
  else{ 
   min_x <- xlim[1]
   max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    min_y <- floor(min(q2.5[(centers < max_x) & (centers > min_x)]  )) - 100
    max_y <- floor(max(q97.5[(centers < max_x) & (centers > min_x)] )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  # artificial max for testing. comment out later
  #max_x <- 90
  
  plot( main = title , xlab = "Age" , ylab = "Hippocampal Volume" ,
        centers , mn,
        ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , 
        xaxt="n" , yaxt="none", type='l')
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 100) , col='grey' , lty=2)
  
  polygon( c( centers , rev(centers) ) , c(q2.5 , rev(q97.5)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  polygon( c( centers , rev(centers) ) , c(q5 , rev(q95)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  polygon( c( centers , rev(centers) ) , c(q10 , rev(q90)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  polygon( c( centers , rev(centers) ) , c(q25 , rev(q75)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  
  lines( centers , q2.5 )   
  lines( centers , q5)
  lines( centers , q10)
  lines( centers , q25)
  lines( centers , mn)
  lines( centers , q75)
  lines( centers , q90)
  lines( centers , q95)
  lines( centers , q97.5)
}


PLOT_NOMOGRAM_COMPARE <- function( bins_1 , bins_2 , hem="bilateral" , col_1="blue" , col_2="red" , title = " " , xlim=NA , ylim=NA , shade=TRUE){
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
  percentiles <- paste("q",LEVELS , suf , sep = "")

  
  if( is.na(xlim[1])){
    min_x <- floor( min( head(bins_1$mid_age,1) , head(bins_2$mid_age,1) )) - 2
    max_x <- floor( max( tail(bins_1$mid_age,1) , tail(bins_2$mid_age,1) )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    # we want the y scale to be outside the most extreme values between the bins being plotted
    # we want the y scale to be calculated within the x limit given above
    bins_1_in_range <- bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , ]
    bins_2_in_range <- bins_2[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , ]
    min_y <- floor(min( c(bins_1_in_range[ , percentiles[1] ] , bins_2_in_range[ , percentiles[1] ]) )) - 100
    max_y <- floor(max( c(bins_1_in_range[ , percentiles[length(LEVELS)] ] , bins_2_in_range[ , percentiles[length(LEVELS)] ] ) )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  
  plot( main = title , bins_1$mid_age , bins_1[, percentiles[floor(length(LEVELS)/2)] ], xlab = "Age" , ylab = "Hippocampal Volume" ,
        type='l' , col=col_1 , lwd=2 , ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , xaxt="n" , yaxt="none")
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 100) , col='grey' , lty=2)
  
  for ( i in percentiles ){
    
    lines(bins_1$mid_age , bins_1[,i] , col=col_1 , lwd=2 ) 
    lines(bins_2$mid_age , bins_2[,i] , col=col_2 , lwd=2 ) 
    if(shade){
    polygon( c( bins_1$mid_age  ,   rev(bins_2$mid_age)  ), 
             c(bins_1[,i], rev(bins_2[,i]) ), 
             col=rgb(0,0,0,0.5) , border = NA , density = 40 )
    }
  }
  
}


PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL <- function( bins_1 , bins_2 , bins_3 , hem="bilateral" , col_1 , col_2 , col_3, title, xlim=NA , ylim=NA , shade=TRUE){
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
  percentiles <- paste("q", LEVELS , suf , sep = "")
  
    if( is.na(xlim[1])){
    min_x <- floor( min( head(bins_1$mid_age,1) , head(bins_2$mid_age,1) , head(bins_3$mid_age,1) )) - 2
    max_x <- floor( max( tail(bins_1$mid_age,1) , tail(bins_2$mid_age,1) , tail(bins_3$mid_age,1) )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    bins_1_in_range <- bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , ]
    bins_2_in_range <- bins_2[ (bins_2$mid_age < max_x) & (bins_2$mid_age > min_x) , ]
    bins_3_in_range <- bins_3[ (bins_3$mid_age < max_x) & (bins_3$mid_age > min_x) , ]
    min_y <- floor(min( c(bins_1_in_range[ , percentiles[1] ] ,
                          bins_2_in_range[ , percentiles[1] ] ,
                          bins_3_in_range[ , percentiles[1] ]) )) - 100
    max_y <- floor(max( c(bins_1_in_range[ , percentiles[length(percentiles)] ] ,
                          bins_2_in_range[ , percentiles[length(percentiles)] ] ,
                          bins_3_in_range[ , percentiles[length(percentiles)] ]) )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  
  plot( main = title , bins_1$mid_age , bins_1[, paste("q50",suf,sep="") ], xlab = "Age" , ylab = "HV" ,
        type='l' , col=col_1 , lwd=2 , ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , xaxt="n" , yaxt="none")

  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 100) , col='grey' , lty=2)
  
  for ( i in percentiles ){
    
    lines(bins_1$mid_age , bins_1[,i] , col=col_1 , lwd=2 ) 
    lines(bins_2$mid_age , bins_2[,i] , col=col_2 , lwd=2 ) 
    lines(bins_3$mid_age , bins_3[,i] , col=col_3 , lwd=2 ) 
    if(shade){
      polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
               c(bins_1[,i], rev(bins_2[,i]) ), 
               col=rgb(0,0,0,0.5) , border = NA , density = 40 )
      
      polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
               c(bins_1[,i], rev(bins_3[,i]) ), 
               col=rgb(0,0,0,0.5) , border = NA , density = 40 )
    }
  }
  
}

PLOT_NOMOGRAM_ADNI <- function( bins , adni , hem =  "bilateral" , title = "" , xlim=NA , ylim=NA ){
  
  adni_AD <- filter(adni , VISCODE.y=="scmri" & DX=="Dementia" )
  adni_MCI <- filter(adni , VISCODE.y=="scmri" & DX=="MCI" )
  adni_CN <- filter(adni , VISCODE.y=="scmri" & DX=="CN" ) 
  adni_NA <- filter(adni , VISCODE.y=="scmri" & is.na(DX) )
  hv_col <- paste( "clean_hv_",hem , sep="")
    
  PLOT_NOMOGRAM( bins , hem = hem , title = title , xlim = xlim , ylim = ylim)
  
  points((adni_AD$TRUE.AGE) , (adni_AD[,hv_col]) , pch=17 , col="red")
  points((adni_MCI$TRUE.AGE) , (adni_MCI[,hv_col]) , pch=16 , col="green")
  points((adni_CN$TRUE.AGE) , (adni_CN[,hv_col]) , pch=15 , col="blue")
  points((adni_NA$TRUE.AGE) , (adni_NA[,hv_col]) , pch=15 , col="grey")
  
}

PLOT_NOMOGRAM_ADNI_LONGITUDINAL <- function(bins , adni , hem="bilateral" , title="" , xlim = NA , ylim=NA ){
  
  PLOT_NOMOGRAM( bins , hem=hem , title=title, xlim = xlim, ylim=ylim)
  table <- adni[ c("RID.x" , "TRUE.AGE" , "clean_hv_left", "clean_hv_right" , "clean_hv_bilateral" , "DX")]
  hv_col <- paste("clean_hv_", hem, sep="")
  IDS <- unique(table$RID.x)
  for (ID in IDS ){
    p_table <- table[ table$RID.x==ID , ]
    p_table <- p_table[ order(p_table$TRUE.AGE) , ]
    if( length(p_table[ , hv_col]) > 1)
      #if(p_table[1,"TRUE.AGE"] < 73 )
      if( sum(!is.na(unique(p_table$DX))) > 1 ){
        lines( p_table$TRUE.AGE , p_table[,hv_col] , col=rgb(0,0,0, alpha = 0.5))
        points((p_table[ p_table$DX=="Dementia","TRUE.AGE"]) , (p_table[ p_table$DX=="Dementia",hv_col]) , pch=17 , col="red")
        points((p_table[ p_table$DX=="MCI","TRUE.AGE"]) , (p_table[ p_table$DX=="MCI",hv_col]) , pch=16 , col="green")
        points((p_table[ p_table$DX=="CN","TRUE.AGE"]) , (p_table[ p_table$DX=="CN",hv_col]), pch=15 , col="blue")
        points((p_table[ is.na(p_table$DX),"TRUE.AGE"]) , (p_table[ is.na(p_table$DX),hv_col]) , pch=15 , col="grey")
        text( x = p_table[1,"TRUE.AGE"] , y = p_table[1,hv_col] , labels = p_table[1,"RID.x"] )
      }
  }
  
}

AVERAGE_PERCENTILE_INTERPOLATE <- function( ADNI , bins , hem  ){
  
  pers <- ( 1:nrow(ADNI) ) * NA
  
  for( i in 1:nrow(ADNI)){
    age <- ADNI$TRUE.AGE[i]
    hv <- ADNI[ i , hv_col ]
    if( !is.na(age) & !is.na(hv)){
      if( (age > min(bins$mid_age)) & (age < max(bins$mid_age)) ){
        
        min_ind <- which.min( abs(age - bins$mid_age) )
        percent.borders <- c(0, LEVELS, 100)
        
        suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
        percentiles <- paste("q",LEVELS , suf , sep = "")
        hv.borders <- bins[ min_ind , percentiles]

        percentile <- 0
        pers[i] <- 0
        if (hv < min(hv.borders)){
          percentile <- 0
          pers[i] <- 0
        } else {
          if (hv > max(hv.borders)){
            percentile <- 100
            pers[i] <- 100
          } else {
            over_hv  <- as.double( hv.borders[ which.min(hv - hv.borders > 0) ] )
            under_hv <- as.double( hv.borders[ which.min(hv - hv.borders > 0) -1 ] )
            over_perc <-  as.double( percent.borders[ which.min(hv - hv.borders > 0) ] )
            under_perc <- as.double( percent.borders[ which.min(hv - hv.borders > 0) -1 ])
            percentile <- (( (hv - under_hv) / (over_hv - under_hv) ) * (over_perc - under_perc) ) + under_perc
            pers[i] <- unlist( (( (hv - under_hv) / (over_hv - under_hv) ) * (over_perc - under_perc) ) + under_perc )
          }
        }
      }
    }
  }
  
  return(pers)
} 

AVERAGE_PERCENTILE <- function( ADNI , bins , hv_col_name , hem = "left" , xlim=NA){
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
  
  if( is.na(xlim[1]) ){
  out_of_range <- ADNI$TRUE.AGE < min(bins$mid_age) | ADNI$TRUE.AGE > max(bins$mid_age)
  }
  else{
    out_of_range <- ADNI$TRUE.AGE < xlim[1] | ADNI$TRUE.AGE > xlim[2]
  }
  indices <- unlist( apply ( as.matrix(ADNI$TRUE.AGE), 1 , function(x) which.min( abs(x - bins$mid_age) ) ) )
  pers <- pnorm( ADNI[ , hv_col_name ] , mean = bins[indices , paste("mean",suf,sep="")] ,
                                         sd = bins[indices , paste("std",suf,sep="")] )

  pers[out_of_range] <- NA 
  return(pers*100)

} 

# Given a table, sort by given column, and return the top and bottom given percent of rows in two separate objects
SEPERATE <- function( Table , Column , Percent){
  
  Table <- Table[ order(Table[,Column]), ]

  Thresh <- nrow(Table) * Percent

  Table_low <- head( Table , Thresh )
  Table_high <- tail( Table , Thresh )
  
  return(list(Table_low = Table_low, Table_high = Table_high))
}

MCI_Analysis <- function( ADNI_table = ADNI_male , hem = "bilateral" , bins = bins_male_gpr , title=""){
  
  #PLOT_NOMOGRAM( bins , hem=hem, xlim = c(30,100), ylim=c(2000,5200) , title = title)
  table <- ADNI_table[ c("RID.x" , "TRUE.AGE" , "clean_hv_left", "clean_hv_right" , "clean_hv_bilateral" , "DX" , "DX.bl")]
  hv_col <- paste("clean_hv_",hem , sep = "")
  IDS <- unique(table$RID.x)
  MCI_table <- data.frame(matrix(ncol = 5, nrow = 0))
  names(MCI_table) <- c("RID" , "age" , "hv" , "percent" , "convert")
  MCI_to_AD <- c()
  MCI_to_MCI <- c()
  for (ID in IDS ){
    p_table <- table[table$RID.x==ID , ]
    p_table <- p_table[ order(p_table$TRUE.AGE) , ]
    if( length(p_table[ , hv_col]) > 1){
      #if( sum(!is.na(unique(p_table$DX))) > 1 )
      if( "MCI" %in% p_table$DX ){
        first_MCI_ind <- match( "MCI" , p_table$DX )
        first_MCI_hv <- p_table[ first_MCI_ind , hv_col]
        first_MCI_perc <- AVERAGE_PERCENTILE( ADNI = p_table[ first_MCI_ind , ] , bins = bins , hem = hem , hv_col_name = hv_col)
        if(p_table[first_MCI_ind,"DX.bl"] == "EMCI")
          if(p_table[first_MCI_ind,"TRUE.AGE"] < 82 ){
            if( "Dementia" %in% p_table$DX ){
              #points( p_table[first_MCI_ind,"TRUE.AGE"] , first_MCI_hv , pch=15 , col="red")
              #text( x = p_table[first_MCI_ind,"TRUE.AGE"] , y = first_MCI_hv+100 , labels = floor(first_MCI_perc) )
              
              MCI_to_AD <- c( MCI_to_AD , first_MCI_perc)
              MCI_table[nrow(MCI_table) + 1,] <- c(unlist(p_table[first_MCI_ind, c("RID.x","TRUE.AGE") ]) ,first_MCI_hv, first_MCI_perc, 1)
            }
            else{
              #points( p_table[first_MCI_ind,"TRUE.AGE"] , first_MCI_hv , pch=15 , col="blue")
              #text( x = p_table[first_MCI_ind,"TRUE.AGE"] , y = first_MCI_hv+100 , labels = floor(first_MCI_perc) )
              
              MCI_to_MCI <- c( MCI_to_MCI , first_MCI_perc)
              MCI_table[nrow(MCI_table) + 1,] <- c(unlist(p_table[first_MCI_ind, c("RID.x","TRUE.AGE") ]) ,first_MCI_hv, first_MCI_perc , 0)
              
              
            }
          }
      }
      
    }
  }
  print( paste("MCI to MCI mean" , mean( MCI_to_MCI , na.rm = TRUE) , "and " , "MCI to AD mean " , mean(MCI_to_AD , na.rm=TRUE) ) )
  return(MCI_table)
}

