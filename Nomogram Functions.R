# File Containing all functions that deal with Nomograms
# By: Mohammed Janahi
#
#


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
#    l.q2.5 : 
#    l.q5 : 
#    ...
#    l.q97.5 :  
#    max_left : maximum left hippocampal volume of samples included in this bin/window
#    mean_right : same as before, but for right HV
#    ...
#    max_right : same as before, but for right HV

# ASSUMPTIONS:
#
WINDOW_ANALYSIS <- function( ukb , percent_samples_per_bin=10 , kernel_width = 20,
                             age_column="AGE_Latest" , 
                             hv_left_column="clean_hv_left", 
                             hv_right_column="clean_hv_right"){
  
  # sort samples by age
  ukb <- ukb[ order(ukb[,age_column]), ]
  
  NUM_SAMPLES <- nrow(ukb)
  
  # calculate the size of each bin
  samples_per_bin <- NUM_SAMPLES / percent_samples_per_bin
  # 1% of samples will determine the start  index of each bin
  one_percent_of_samples <- NUM_SAMPLES/100
  
  # initialize the output table
  bins <- data.frame(min_age=double(), max_age=double(), mid_age=double(), 
                     
                     mean_left=double(), std_left=double(),
                     min_left=double(), 
                     l.q2.5=double(), l.q5=double(), l.q10=double(), l.q25=double(), l.q50=double(), l.q75=double(), l.q90=double(), l.q95=double(), l.q97.5=double(), 
                     max_left=double(),
                     
                     mean_right=double(), std_right=double(),
                     min_right=double(), 
                     r.q2.5=double(), r.q5=double(), r.q10=double(), r.q25=double(), r.q50=double(), r.q75=double(), r.q90=double(), r.q95=double(), r.q97.5=double(),
                     max_right=double()
  )
  # loop from beginning, increment by 1% of samples each time, stop before the size of the last bin exceeds num of samples
  for (i in seq( 1 , (NUM_SAMPLES - samples_per_bin) , one_percent_of_samples) ){
    # get samples on this window
    window_of_samples_left  <- ukb[i:(i+samples_per_bin) , c(age_column, hv_left_column)]
    window_of_samples_right <- ukb[i:(i+samples_per_bin) , c(age_column, hv_right_column)]
    
    win_ages      <- window_of_samples_left[ , c(age_column) ]
    win_min_age   <- min(win_ages  , na.rm = TRUE)
    win_max_age   <- max(win_ages  , na.rm = TRUE)
    win_mean_age  <- mean(win_ages , na.rm = TRUE)
    
    win_hv_left   <- window_of_samples_left[,hv_left_column]
    win_hv_right  <- window_of_samples_right[,hv_right_column]
   
    win_mean_hv_left   <- mean(win_hv_left, na.rm = TRUE) 
    win_std_hv_left    <- sd(win_hv_left , na.rm = TRUE)
    win_min_hv_left    <- min(win_hv_left , na.rm = TRUE)
    win_max_hv_left    <- max(win_hv_left , na.rm = TRUE)
    win_quantiles_left <- quantile( (win_hv_left) , probs = c(0.025, 0.05 , 0.1 , 0.25 , 0.5 , 0.75 , 0.9 , 0.95, 0.975) , na.rm = TRUE) 
    
    win_mean_hv_right  <- mean(win_hv_right, na.rm = TRUE) 
    win_std_hv_right   <- sd(win_hv_right , na.rm = TRUE)
    win_min_hv_right   <- min(win_hv_right , na.rm = TRUE)
    win_max_hv_right   <- max(win_hv_right , na.rm = TRUE)
    win_quantiles_right <- quantile( (win_hv_right) , probs = c(0.025, 0.05 , 0.1 , 0.25 , 0.5 , 0.75 , 0.9 , 0.95, 0.975) , na.rm = TRUE) 
    
    row <- unlist(list( win_min_age , win_max_age , win_mean_age ,
                        win_mean_hv_left , win_std_hv_left ,
                        win_min_hv_left ,
                        win_quantiles_left , 
                        win_max_hv_left ,
                        win_mean_hv_right , win_std_hv_right,
                        win_min_hv_right ,
                        win_quantiles_right ,
                        win_max_hv_right  ))
    
    bins[nrow(bins) + 1 , ]  = row
  }
  
  # Gaussian smoothing
  for (col in names(bins)){
    bins[ , col ] <- smth( bins[, col ],  method = 'gaussian', window=kernel_width)
  }
  # return the non-na rows
  return(bins[!is.na(bins$min_age),])
}


GPR_ANALYSIS <- function(ukb , XX=NA , age_column="AGE_Latest" , 
                         hv_left_column="clean_hv_left", hv_right_column="clean_hv_right"){
  if(!is.na(XX[1])){
    L = nrow(XX)
  }
  else{
    L = nrow(ukb)
    XX <- matrix(seq( 30, 100, length = L) , ncol=1)
  }
  
  # initialize the output table
  bins <- data.frame(mid_age=double(L),
                     
                     mean_left=double(L), std_left=double(L),
                     min_left=double(L), 
                     l.q2.5=double(L), l.q5=double(L), l.q10=double(L), l.q25=double(L), l.q50=double(L), 
                     l.q75=double(L), l.q90=double(L), l.q95=double(L), l.q97.5=double(L), 
                     max_left=double(L),
                     
                     mean_right=double(L), std_right=double(L),
                     min_right=double(L), 
                     r.q2.5=double(L), r.q5=double(L), r.q10=double(L), r.q25=double(L), r.q50=double(L), 
                     r.q75=double(L), r.q90=double(L), r.q95=double(L), r.q97.5=double(L),
                     max_right=double(L)
  )
  
  X = ukb[ , age_column]
  y = ukb[ , hv_left_column]
  
  X = X[!is.na(y)]
  y = y[!is.na(y)]
  
  bins$mid_age <- c(XX)
  
  gpi <- newGPsep(matrix(X , ncol=1) ,y , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  p <- predGPsep( gpi , XX )
  
  deleteGPsep(gpi)
  
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
  
  X = ukb[ , age_column]
  y = ukb[ , hv_right_column]
  
  X = X[!is.na(y)]
  y = y[!is.na(y)]
  
  
  gpi <- newGPsep(matrix(X , ncol=1) ,y , d=0.1, g=0.1*var(y) , dK=TRUE) 
  eps <- sqrt(.Machine$double.eps)
  mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  p <- predGPsep( gpi , XX )
  deleteGPsep(gpi)
  
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
  
  return(bins)
}

GPR_ANALYSIS_PARALLEL <- function(ukb , XX=NA , age_column="AGE_Latest" , 
                         hv_left_column="clean_hv_left", hv_right_column="clean_hv_right"){
  if(!is.na(XX[1])){
    L = nrow(XX)
  }
  else{
    L = nrow(ukb)
    XX <- matrix(seq( 30, 100, length = L) , ncol=1)
  }
  
  # initialize the output table
  bins <- data.frame(mid_age=double(L),
                     
                     mean_left=double(L), std_left=double(L),
                     min_left=double(L), 
                     l.q2.5=double(L), l.q5=double(L), l.q10=double(L), l.q25=double(L), l.q50=double(L), 
                     l.q75=double(L), l.q90=double(L), l.q95=double(L), l.q97.5=double(L), 
                     max_left=double(L),
                     
                     mean_right=double(L), std_right=double(L),
                     min_right=double(L), 
                     r.q2.5=double(L), r.q5=double(L), r.q10=double(L), r.q25=double(L), r.q50=double(L), 
                     r.q75=double(L), r.q90=double(L), r.q95=double(L), r.q97.5=double(L),
                     max_right=double(L)
  )
  
  X = ukb[ , age_column]
  y = ukb[ , hv_left_column]
  
  X = X[!is.na(y)]
  y = y[!is.na(y)]
  
  bins$mid_age <- c(XX)
  
  p <- aGP( X = matrix(X , ncol=1) , Z = y , XX = XX, d=0.1, g=0.1*var(y) , omp.threads = 16 , method="mspe" , end = 30)
  
  
  bins$mean_left <- p$mean
  bins$std_left <- sqrt(p$var)
  bins$l.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(p$var))
  bins$l.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(p$var))
  bins$l.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(p$var))
  bins$min_left <- -Inf
  bins$max_left <- Inf
  
  X = ukb[ , age_column]
  y = ukb[ , hv_right_column]
  
  X = X[!is.na(y)]
  y = y[!is.na(y)]
  
  
  p <- aGP( X = matrix(X , ncol=1) , Z = y , XX = XX, d=0.1, g=0.1*var(y) , omp.threads = 16)
  
  
  bins$mean_right <- p$mean
  bins$std_right <- sqrt(p$var)
  bins$r.q2.5  <- qnorm(0.025, mean = p$mean, sd = sqrt(p$var))
  bins$r.q5    <- qnorm(0.05,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q10   <- qnorm(0.10,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q25   <- qnorm(0.25,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q50   <- qnorm(0.50,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q75   <- qnorm(0.75,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q90   <- qnorm(0.90,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q95   <- qnorm(0.95,  mean = p$mean, sd = sqrt(p$var))
  bins$r.q97.5 <- qnorm(0.975, mean = p$mean, sd = sqrt(p$var))
  bins$min_right <- -Inf
  bins$max_right <- Inf
  
  return(bins)
}

GPR_ANALYSIS_TEST <- function(ukb , XX=NA , age_column="AGE_Latest" , 
                         hv_left_column="clean_hv_left", hv_right_column="clean_hv_right"){
  L = nrow(ukb)
  if(is.na(XX[1])){
    XX_male <- c( seq( 30, 100, length = L) ,  (seq( 30, 100, length = L)*0 + 1) )
    dim(XX_male) <- c(L,2)
    XX_female <- c( seq( 30, 100, length = L) ,  (seq( 30, 100, length = L)*0) )
    dim(XX_female) <- c(L,2)
  }
  
  # initialize the output table
  bins_male <- data.frame(mid_age=double(L),
                          
                          mean_left=double(L), std_left=double(L),
                          min_left=double(L), 
                          l.q2.5=double(L), l.q5=double(L), l.q10=double(L), l.q25=double(L), l.q50=double(L), 
                          l.q75=double(L), l.q90=double(L), l.q95=double(L), l.q97.5=double(L), 
                          max_left=double(L),
                          
                          mean_right=double(L), std_right=double(L),
                          min_right=double(L), 
                          r.q2.5=double(L), r.q5=double(L), r.q10=double(L), r.q25=double(L), r.q50=double(L), 
                          r.q75=double(L), r.q90=double(L), r.q95=double(L), r.q97.5=double(L),
                          max_right=double(L)
  )
  
  
  # initialize the output table
  bins_female <- data.frame(mid_age=double(L),
                          
                          mean_left=double(L), std_left=double(L),
                          min_left=double(L), 
                          l.q2.5=double(L), l.q5=double(L), l.q10=double(L), l.q25=double(L), l.q50=double(L), 
                          l.q75=double(L), l.q90=double(L), l.q95=double(L), l.q97.5=double(L), 
                          max_left=double(L),
                          
                          mean_right=double(L), std_right=double(L),
                          min_right=double(L), 
                          r.q2.5=double(L), r.q5=double(L), r.q10=double(L), r.q25=double(L), r.q50=double(L), 
                          r.q75=double(L), r.q90=double(L), r.q95=double(L), r.q97.5=double(L),
                          max_right=double(L)
  )
  
  X = ukb[ , c(age_column,"SEX")]
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
  
  X = ukb[ , c(age_column,"SEX")]
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
  
  return( list(bins_male = bins_male,bins_female = bins_female) )
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


#Assumptions:
# - Bins have the exact same mid_age column entries (or very close)
# - Bins are constructed as defined in the ANALYSIS functions (they have columns l.q* and r.q* )
NOMOGRAM_DIFF <- function( bins_1 , bins_2 , hem=NA ){
  error = 0
  
  if(is.na(hem) )
    cols = c("l.q2.5","l.q5","l.q10","l.q25","l.q50","l.q75","l.q90","l.q95","l.q97.5",
             "r.q2.5","r.q5","r.q10","r.q25","r.q50","r.q75","r.q90","r.q95","r.q97.5" )
  else{
    if(hem == 0)
      cols = c("l.q2.5","l.q5","l.q10","l.q25","l.q50","l.q75","l.q90","l.q95","l.q97.5")
    if(hem == 1)
      cols = c("r.q2.5","r.q5","r.q10","r.q25","r.q50","r.q75","r.q90","r.q95","r.q97.5" )
  }
  
  for ( col_name in cols )
    error = error + sum( abs( bins_1[ , col_name] - bins_2[ , col_name] ) , na.rm = TRUE )
  
  return(error)
}

NOMOGRAM_DIFF_INTERPOLATE <- function( bins_1 , bins_2 , hem=NA ){
  error = 0
  
    
  if(is.na(hem) )
    cols = c("l.q2.5","l.q5","l.q10","l.q25","l.q50","l.q75","l.q90","l.q95","l.q97.5",
             "r.q2.5","r.q5","r.q10","r.q25","r.q50","r.q75","r.q90","r.q95","r.q97.5" )
  else{
    if(hem == 0)
      cols = c("l.q2.5","l.q5","l.q10","l.q25","l.q50","l.q75","l.q90","l.q95","l.q97.5")
    if(hem == 1)
      cols = c("r.q2.5","r.q5","r.q10","r.q25","r.q50","r.q75","r.q90","r.q95","r.q97.5" )
  }
  
  age_range_min <- max(min(bins_1$mid_age),min(bins_2$mid_age))
  age_range_max <- min(max(bins_1$mid_age),max(bins_2$mid_age))
  
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
    #error = error + sum( abs( bins_1[ , col_name] - bins_2[ , col_name] ) , na.rm = TRUE )
  max_age_range <- max(max(bins_1$mid_age),max(bins_2$mid_age)) - min(min(bins_1$mid_age),min(bins_2$mid_age))
  min_age_range <- (age_range_max - age_range_min)
  overlap_percent <- min_age_range / max_age_range 
  print(paste("AGE RANGE OVERLAP:" , overlap_percent))
  return(error)
}

# Function to plot nomograms given a bins object
#
PLOT_NOMOGRAM <- function( bins , hem = 0 , title = " " , xlim = NA , ylim = NA ){
  
  centers <- bins$mid_age
  
  if( hem == 0 ){
    mn <- bins$mean_left
    q2.5 <- bins$l.q2.5
    q5   <- bins$l.q5    
    q10  <- bins$l.q10 
    q25  <- bins$l.q25 
    q50  <- bins$l.q50 
    q75  <- bins$l.q75 
    q90  <- bins$l.q90 
    q95  <- bins$l.q95 
    q97.5<- bins$l.q97.5 
  }
  if( hem == 1 ){
    mn <- bins$mean_right 
    q2.5 <- bins$r.q2.5
    q5   <- bins$r.q5 
    q10  <- bins$r.q10 
    q25  <- bins$r.q25 
    q50  <- bins$r.q50 
    q75  <- bins$r.q75
    q90  <- bins$r.q90 
    q95  <- bins$r.q95 
    q97.5<- bins$r.q97.5 
  }
  if( hem == 2 ){
    mn <- bins$mean_left  + bins$mean_right 
    q2.5 <- bins$l.q2.5  + bins$r.q2.5
    q5   <- bins$l.q5    + bins$r.q5 
    q10  <- bins$l.q10   + bins$r.q10 
    q25  <- bins$l.q25   + bins$r.q25
    q50  <- bins$l.q50   + bins$r.q50 
    q75  <- bins$l.q75   + bins$r.q75 
    q90  <- bins$l.q90   + bins$r.q90
    q95  <- bins$l.q95   + bins$r.q95 
    q97.5<- bins$l.q97.5 + bins$r.q97.5
  }

  

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
  
  plot( main = title , xlab = "Age" , ylab = "HV" ,
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


PLOT_NOMOGRAM_COMPARE <- function( bins_1 , bins_2 , left=0 , col_1="blue" , col_2="red" , title = " " , xlim=NA , ylim=NA ){
  
  if (left==1) percentiles <- c("l.q97.5","l.q95","l.q90","l.q75","l.q50","l.q25","l.q10","l.q5","l.q2.5")
  if (left==0) percentiles <- c("r.q97.5","r.q95","r.q90","r.q75","r.q50","r.q25","r.q10","r.q5","r.q2.5")
  
  
  if( is.na(xlim[1])){
    min_x <- floor( min( head(bins_1$mid_age,1) , head(bins_2$mid_age,1) )) - 2
    max_x <- floor( max( tail(bins_1$mid_age,1) , tail(bins_2$mid_age,1) )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    # we want the y scale to be outside the most extrme values between the bins being plotted
    # we want the y scale to be calculated within the x limit given above
    min_y <- floor(min( c(bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , percentiles[9] ] , 
                          bins_2[ (bins_2$mid_age < max_x) & (bins_2$mid_age > min_x) , percentiles[9] ]) )) - 100
    max_y <- floor(max( c(bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , percentiles[1] ] ,
                          bins_2[ (bins_2$mid_age < max_x) & (bins_2$mid_age > min_x) , percentiles[1] ] ) )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  
  plot( main = title , bins_1$mid_age , bins_1[, percentiles[5] ], xlab = "Age" , ylab = "HV" ,
        type='l' , col=col_1 , lwd=2 , ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , xaxt="n" , yaxt="none")
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 100) , col='grey' , lty=2)
  
  for ( i in percentiles ){
    
    lines(bins_1$mid_age , bins_1[,i] , col=col_1 , lwd=2 ) 
    lines(bins_2$mid_age , bins_2[,i] , col=col_2 , lwd=2 ) 
    
    polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
             c(bins_1[,i], rev(bins_2[,i]) ), 
             col=rgb(0,0,0,0.5) , border = NA , density = 40 )
  }
  
}


PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL <- function( bins_1 , bins_2 , bins_3 , left , col_1 , col_2 , col_3, title, xlim=NA , ylim=NA){
  
  if (left==1) percentiles <- c("l.q97.5","l.q95","l.q90","l.q75","l.q50","l.q25","l.q10","l.q5","l.q2.5")
  if (left==0) percentiles <- c("r.q97.5","r.q95","r.q90","r.q75","r.q50","r.q25","r.q10","r.q5","r.q2.5")

    if( is.na(xlim[1])){
    min_x <- floor( min( head(bins_1$mid_age,1) , head(bins_2$mid_age,1) , head(bins_3$mid_age,1) )) - 2
    max_x <- floor( max( tail(bins_1$mid_age,1) , tail(bins_2$mid_age,1) , tail(bins_3$mid_age,1) )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    min_y <- floor(min( c(bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , percentiles[9] ] ,
                          bins_2[ (bins_2$mid_age < max_x) & (bins_2$mid_age > min_x) , percentiles[9] ] ,
                          bins_3[ (bins_3$mid_age < max_x) & (bins_3$mid_age > min_x), percentiles[9] ]) )) - 100
    max_y <- floor(max( c(bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , percentiles[1] ] ,
                          bins_2[ (bins_2$mid_age < max_x) & (bins_2$mid_age > min_x) , percentiles[1] ] ,
                          bins_3[ (bins_3$mid_age < max_x) & (bins_3$mid_age > min_x), percentiles[1] ] ) )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  
  plot( main = title , bins_1$mid_age , bins_1[, percentiles[5] ], xlab = "Age" , ylab = "HV" ,
        type='l' , col=col_1 , lwd=2 , ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , xaxt="n" , yaxt="none")

  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 100) , col='grey' , lty=2)
  
  for ( i in percentiles ){
    
    lines(bins_1$mid_age , bins_1[,i] , col=col_1 , lwd=2 ) 
    lines(bins_2$mid_age , bins_2[,i] , col=col_2 , lwd=2 ) 
    lines(bins_3$mid_age , bins_3[,i] , col=col_3 , lwd=2 ) 
    
    polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
             c(bins_1[,i], rev(bins_2[,i]) ), 
             col=rgb(0,0,0,0.5) , border = NA , density = 40 )
    
    polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
             c(bins_1[,i], rev(bins_3[,i]) ), 
             col=rgb(0,0,0,0.5) , border = NA , density = 40 )
    
  }
  
}


AVERAGE_PERCENTILE_INTERPOLATE <- function( ADNI , bins , hv_col_name , left ){
  
  pers <- ( 1:nrow(ADNI) ) * NA
  
  for( i in 1:nrow(ADNI)){
    age <- ADNI$TRUE.AGE[i]
    hv <- ADNI[ i , hv_col_name ]
    if( !is.na(age) & !is.na(hv)){
      if( (age > min(bins$mid_age)) & (age < max(bins$mid_age)) ){
        min_ind <- which.min( abs(age - bins$mid_age) )
        percent.borders <- c(0, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 100)
        if(left == 1) hv.borders <- bins[min_ind , 19:29] # indices of columns for left hv (better to mention by column name not index for modularity)
        if(left == 0) hv.borders <- bins[min_ind , 32:42]
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
  #return( total_perc / sum( !is.na( ADNI$AGE ) ) )
} 

AVERAGE_PERCENTILE <- function( ADNI , bins , hv_col_name , left , xlim=NA){
  if( is.na(xlim[1]) ){
  out_of_range <- ADNI$TRUE.AGE < min(bins$mid_age) | ADNI$TRUE.AGE > max(bins$mid_age)
  }
  else{
    out_of_range <- ADNI$TRUE.AGE < xlim[1] | ADNI$TRUE.AGE > xlim[2]
  }
  indices <- unlist( apply ( as.matrix(ADNI$TRUE.AGE), 1 , function(x) which.min( abs(x - bins$mid_age) ) ) )
  
  if(left == 1) pers <- pnorm( ADNI[ , hv_col_name ] , mean = bins[indices , "mean_left"] , sd = bins[indices , "std_left"] )
  if(left == 0) pers <- pnorm( ADNI[ , hv_col_name ] , mean = bins[indices , "mean_right"] , sd = bins[indices , "std_right"] )
  
  pers[out_of_range] <- NA 
  return(pers*100)

} 

# Given a table, sort by given column, and return the top and bottom given percent of rows in two separate objects
SEPERATE <- function( Table , Column , Percent){
  
  Table <- Table[ order(Table[,Column]), ]

  Thresh_male <- nrow(Table) * Percent

  Table_low <- head( Table , Thresh_male )
  Table_high <- tail( Table , Thresh_male )
  
  return(list(Table_low = Table_low, Table_high = Table_high))
}


THRESHOLD_EXPLORE <- function( ukb_img_ex_outliers_male , ukb_img_ex_outliers_female ){ 
  # this is exploration of where is the best threshold/ percentile to cut samples at 
  percentiles_l <- c("l.q97.5","l.q95","l.q90","l.q75","l.q50","l.q25","l.q10","l.q5","l.q2.5")
  percentiles_r <- c("r.q97.5","r.q95","r.q90","r.q75","r.q50","r.q25","r.q10","r.q5","r.q2.5")
  thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.06','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                  'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')
  colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
              '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  
  #plot(0 , xlim = c(0,length(50:10)) , ylim = c(0,900) , col = "#ffffff" , xaxt = "n" , xlab="Extremes Percent" , ylab = "MAE")
  #axis(1, at=0:40,labels=50:10)
  par(mfrow=c(4,4))
  
  errors <- matrix(0, length(thresholds) , length(50:10))
  for( perc in seq(50,10,-4)){
    ind <- 1
    for (PRS_Thresh in c("PRS_TH_1e.08" , "PRS_TH_1e.06" , "PRS_TH_0.2" , "PRS_TH_1") ){
      
      ukb_img_ex_outliers_male <- ukb_img_ex_outliers_male[ order(ukb_img_ex_outliers_male[,PRS_Thresh]), ]
      ukb_img_ex_outliers_female <- ukb_img_ex_outliers_female[ order(ukb_img_ex_outliers_female[,PRS_Thresh]), ]
      
      Thresh_male <- nrow(ukb_img_ex_outliers_male) * (perc/100)
      Thresh_female <- nrow(ukb_img_ex_outliers_female) * (perc/100)
      
      ukb_img_ex_outliers_male_lower <- head( ukb_img_ex_outliers_male , Thresh_male )
      ukb_img_ex_outliers_male_upper <- tail( ukb_img_ex_outliers_male , Thresh_male )
      
      ukb_img_ex_outliers_female_lower <- head( ukb_img_ex_outliers_female , Thresh_female )
      ukb_img_ex_outliers_female_upper <- tail( ukb_img_ex_outliers_female , Thresh_female )
      
      bins_male      <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male )
      bins_male_low  <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male_lower )
      bins_male_high <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male_upper ) 
      
      bins_female      <- WINDOW_ANALYSIS( ukb_img_ex_outliers_female )
      bins_female_low  <- WINDOW_ANALYSIS( ukb_img_ex_outliers_female_lower )
      bins_female_high <- WINDOW_ANALYSIS( ukb_img_ex_outliers_female_upper ) 
      
      hist(ukb_img_ex_outliers_male_upper$AGE_Latest , 100 , main = "male UPPER" , xlab = "age" , ylab = paste("P:",perc , "% T:" , substr(PRS_Thresh , 8 , 100) ) )
      hist(ukb_img_ex_outliers_male_lower$AGE_Latest , 100 , main = "male LOWER" , xlab = "age" )
      hist(ukb_img_ex_outliers_female_upper$AGE_Latest , 100 , main = "FEmale UPPER" , xlab = "age" )
      hist(ukb_img_ex_outliers_female_lower$AGE_Latest , 100 , main = "FEmale LOWER" , xlab = "age" )
      
      x <- 0
      for (i in c(percentiles_l,percentiles_r) ){
        x <- x + mean(abs(bins_male[,i] - bins_male_low[,i]) , na.rm = TRUE)
        x <- x + mean(abs(bins_male[,i] - bins_male_high[,i]) , na.rm = TRUE)
        x <- x + mean(abs(bins_female[,i] - bins_female_low[,i]) , na.rm = TRUE)
        x <- x + mean(abs(bins_female[,i] - bins_female_high[,i]) , na.rm = TRUE)
      }
      print(paste( "Percent: " , perc, " Threshold: " , PRS_Thresh, " MAE: " , (x/(length(percentiles)*4) )))
      points(ind , x/(length(percentiles)*4)  , col = colors[match(PRS_Thresh ,thresholds)] , pch = 19)
      errors[ ind , (51 - perc) ] <- x/(length(percentiles)*4)
      ind <- ind + 1
    }
    
  }
}
