# Exploration Functions
#
#

#
# Test PRS built with different SNPs / thresholds
#
EXPLORE_THRESHOLD <- function( ukb_img_ex_outliers_male , ukb_img_ex_outliers_female ){ 
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
      hist(ukb_img_ex_outliers_female_upper$AGE_Latest , 100 , main = "Female UPPER" , xlab = "age" )
      hist(ukb_img_ex_outliers_female_lower$AGE_Latest , 100 , main = "Female LOWER" , xlab = "age" )
      
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



#
# Look at GPR Nomograms built with UKB samples plus the CN samples from ADNI
#
# Results:
#    Adding them improves the older end of the Nomograms
EXPLORE_UKB_ADNI_GPR <- function( ukb_img_ex_outliers_male , ADNI_filter_male_CN ,
                                  ukb_img_ex_outliers_female , ADNI_filter_female_CN,
                                  ukb_img_ex_outliers_male_high , ADNI_filter_male_high_CN,
                                  ukb_img_ex_outliers_female_high , ADNI_filter_female_high_CN,
                                  ukb_img_ex_outliers_male_low, ADNI_filter_male_low_CN,
                                  ukb_img_ex_outliers_female_low, ADNI_filter_female_low_CN ){
  
  data <- rbind( ukb_img_ex_outliers_male[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
                 setNames( ADNI_filter_male_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                           c("AGE_Latest","clean_hv_left","clean_hv_right")) )
  bins_male_gpr_accum <- GPR_ANALYSIS(data)
  
  data <- rbind( ukb_img_ex_outliers_female[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
                 setNames( ADNI_filter_female_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                           c("AGE_Latest","clean_hv_left","clean_hv_right")) )
  bins_female_gpr_accum <- GPR_ANALYSIS(data)
  
  data <- rbind( ukb_img_ex_outliers_male_high[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
                 setNames( ADNI_filter_male_high_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                           c("AGE_Latest","clean_hv_left","clean_hv_right")) )
  bins_male_gpr_accum_high <- GPR_ANALYSIS(data)
  
  data <- rbind( ukb_img_ex_outliers_female_high[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
                 setNames( ADNI_filter_female_high_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                           c("AGE_Latest","clean_hv_left","clean_hv_right")) )
  bins_female_gpr_accum_high <- GPR_ANALYSIS(data)
  
  data <- rbind( ukb_img_ex_outliers_male_low[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
                 setNames( ADNI_filter_male_low_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                           c("AGE_Latest","clean_hv_left","clean_hv_right")) )
  bins_male_gpr_accum_low <- GPR_ANALYSIS(data)
  
  data <- rbind( ukb_img_ex_outliers_female_low[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
                 setNames( ADNI_filter_female_low_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                           c("AGE_Latest","clean_hv_left","clean_hv_right")) )
  bins_female_gpr_accum_low <- GPR_ANALYSIS(data)
  
  
  
  
}



#
# Explore how far we can rely on the nomogram at its ends
#
#
EXPLORE_GPR_LIMITS <- function(ukb_table){
  
  # Test 1: include 1000 samples of people below age 60/70/80/90. How long before it goes out of wack?
#  par(mfrow=c(2,2))
  
#  table_60s <- (ukb_table[(ukb_table$AGE_Latest < 60) , ])
  #  table_60s_samples <- table_60s[ sample(1:nrow(table_60s), 1000, replace = FALSE) , ]
  #  bins_60s <- GPR_ANALYSIS(table_60s_samples)
  #PLOT_NOMOGRAM(bins_60s)
  #abline(v=max(table_60s_samples$AGE_Latest) , col="black")
  #abline(v=min(table_60s_samples$AGE_Latest) , col="black")
  
  #table_70s <- (ukb_table[(ukb_table$AGE_Latest < 70) , ])
  #table_70s_samples <- table_70s[ sample(1:nrow(table_70s), 1000, replace = FALSE) , ]
  #bins_70s <- GPR_ANALYSIS(table_70s_samples)
  #PLOT_NOMOGRAM(bins_70s)
  #abline(v=max(table_70s_samples$AGE_Latest) , col="black")
  #abline(v=min(table_70s_samples$AGE_Latest) , col="black")
  
  #table_80s <- (ukb_table[(ukb_table$AGE_Latest < 80) , ])
  #table_80s_samples <- table_80s[ sample(1:nrow(table_80s), 1000, replace = FALSE) , ]
  #bins_80s <- GPR_ANALYSIS(table_80s_samples)
  #PLOT_NOMOGRAM(bins_80s)
  #abline(v=max(table_80s_samples$AGE_Latest) , col="black")
  #abline(v=min(table_80s_samples$AGE_Latest) , col="black")
  
  #table_90s <- (ukb_table[(ukb_table$AGE_Latest < 90) , ])
  #table_90s_samples <- table_90s[ sample(1:nrow(table_90s), 1000, replace = FALSE) , ]
  #bins_90s <- GPR_ANALYSIS(table_90s_samples)
  #PLOT_NOMOGRAM(bins_90s)
  #abline(v=max(table_90s_samples$AGE_Latest) , col="black")
  #abline(v=min(table_90s_samples$AGE_Latest) , col="black")
  
  # Test 2: take overlapping years, 600 samples between ages, check overlapping areas for how long before it goes out of wack
  indices <- floor(seq(1 , nrow(ukb_male), nrow(ukb_male)/1500 ))
  ukb_male <- ukb_male[ order(ukb_male[,"AGE_Latest"]), ]
  full_table <- ukb_male[indices,]
  
  aged_table <- full_table[ (full_table$AGE_Latest >= 40 & full_table$AGE_Latest <= 60), ]
  table_1 <- aged_table[ sample( 1:nrow(aged_table), 1000 , replace = TRUE), ]
  PLOT_NOMOGRAM( bins_1 <- GPR_ANALYSIS(table_1) )
  abline(v=max(table_1$AGE_Latest) , col="black")
  abline(v=min(table_1$AGE_Latest) , col="black")
  
  aged_table <- full_table[ (full_table$AGE_Latest >= 50 & full_table$AGE_Latest <= 70), ]
  table_2 <- aged_table[ sample( 1:nrow(aged_table), 1000 , replace = TRUE), ]
  PLOT_NOMOGRAM( bins_2 <- GPR_ANALYSIS(table_2) )
  abline(v=max(table_2$AGE_Latest) , col="black")
  abline(v=min(table_2$AGE_Latest) , col="black")
  
  aged_table <- full_table[ (full_table$AGE_Latest >= 60 & full_table$AGE_Latest <= 80), ]
  table_3 <- aged_table[ sample( 1:nrow(aged_table), 1000 , replace = TRUE), ]
  PLOT_NOMOGRAM( bins_3 <- GPR_ANALYSIS(table_3) )
  abline(v=max(table_3$AGE_Latest) , col="black")
  abline(v=min(table_3$AGE_Latest) , col="black")
  
  PLOT_NOMOGRAM_COMPARE(bins_1 , bins_2 , shade = FALSE , ylim = c(3000,6000))
  abline(v=max(table_2$AGE_Latest) , col="red")
  abline(v=max(table_1$AGE_Latest) , col="blue")
  
  PLOT_NOMOGRAM_COMPARE(bins_2 , bins_3 , shade = FALSE , ylim = c(3000,6000))
  abline(v=max(table_2$AGE_Latest) , col="blue")
  abline(v=max(table_3$AGE_Latest) , col="red")
  
  PLOT_NOMOGRAM_COMPARE(bins_1 , bins_3 , shade = FALSE , ylim = c(3000,6000))
  abline(v=max(table_3$AGE_Latest) , col="blue")
  abline(v=max(table_1$AGE_Latest) , col="red")
  
  
  # Test 3: does sample size make a difference to reliability of GPR beyond age range? 
  # keep max/min age range constant and try 600/1000/2000/5000 samples
  ukb_table <- ukb_table[ order(ukb_table[,"AGE_Latest"]), ]
  table_all <- (ukb_table[(ukb_table$AGE_Latest < 70) & (ukb_table$AGE_Latest > 60) , ])
  min_aged <- head(table_all , n = 1)
  max_aged <- tail(table_all , n = 1)
  
  table_600_samples <- rbind( min_aged , table_all[ sample(2:(nrow(table_all)-1), 598, replace = FALSE) , ] , max_aged)
  bins_600 <- GPR_ANALYSIS(table_600_samples)
  
  table_1000_samples <- rbind( min_aged , table_all[ sample(2:(nrow(table_all)-1), 998, replace = FALSE) , ] , max_aged)
  bins_1000 <- GPR_ANALYSIS(table_1000_samples)
  
  table_2000_samples <- rbind( min_aged , table_all[ sample(2:(nrow(table_all)-1), 1998, replace = FALSE) , ] , max_aged)
  bins_2000 <- GPR_ANALYSIS(table_2000_samples)
  
  table_5000_samples <- rbind( min_aged , table_all[ sample(2:(nrow(table_all)-1), 4998, replace = FALSE) , ] , max_aged)
  bins_5000 <- GPR_ANALYSIS(table_5000_samples)
  
  
  PLOT_NOMOGRAM(bins_600)
  abline(v=max(table_600_samples$AGE_Latest) , col="blue")
  abline(v=min(table_600_samples$AGE_Latest) , col="blue")
  
  PLOT_NOMOGRAM(bins_1000)
  abline(v=max(table_1000_samples$AGE_Latest) , col="blue")
  abline(v=min(table_1000_samples$AGE_Latest) , col="blue")
  
  PLOT_NOMOGRAM(bins_2000)
  abline(v=max(table_2000_samples$AGE_Latest) , col="blue")
  abline(v=min(table_2000_samples$AGE_Latest) , col="blue")
  
  PLOT_NOMOGRAM(bins_5000)
  abline(v=max(table_5000_samples$AGE_Latest) , col="blue")
  abline(v=min(table_5000_samples$AGE_Latest) , col="blue")
  
}


#
# Explore combining L/R and M/F samples and putting them all in one GPR
#
#
EXPLORE_COMBINED_GPR <- function(ukb_img_ex_outliers_male ){
  
  # exploring putting multiple dimensions into the GPR model
  ukb_img_ex_outliers_male <- ukb_img_ex_outliers_male[ order(ukb_img_ex_outliers_male[,"PRS_TH_1"]), ]
  third <- nrow(ukb_img_ex_outliers_male) * 0.3
  ukb_img_ex_outliers_male$PRS_TH_1_Levels = 0
  ukb_img_ex_outliers_male[ 1:third , "PRS_TH_1_Levels"] = 1
  ukb_img_ex_outliers_male[ nrow(ukb_img_ex_outliers_male):(nrow(ukb_img_ex_outliers_male)-third) , "PRS_TH_1_Levels"] = -1
  
  samples_high <- (ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$PRS_TH_1_Levels == 1 , ])[1:900,]
  samples_low <- (ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$PRS_TH_1_Levels == -1 , ])[1:900,]
  
  samples_combined <- rbind( (ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$PRS_TH_1_Levels == 1 , ])[1:900,] , 
                             (ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$PRS_TH_1_Levels == 0 , ])[1:900,] ,
                             (ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$PRS_TH_1_Levels == -1 , ])[1:900,] )
  
  bins_high <- GPR_ANALYSIS(samples_high)
  bins_low <- GPR_ANALYSIS(samples_low)
  bins_x <- GPR_ANALYSIS_TEST(samples_combined)
  bins_high_combined <- bins_x$bins_high
  bins_low_combined <- bins_x$bins_low
  
  par(mfrow=c(2,2))
  PLOT_NOMOGRAM(bins_high)
  PLOT_NOMOGRAM(bins_low)
  PLOT_NOMOGRAM(bins_high_combined)
  PLOT_NOMOGRAM(bins_low_combined)
  
  # Explore combining L/R hv into one column
  ukb_male_1000 <- ukb_img_ex_outliers_male[1:1000 , ]
  ukb_male_1000_LR <- ukb_male_1000
  ukb_male_1000_R <- ukb_male_1000
  
  ukb_male_1000_LR[ , "L/R" ] <- 1
  ukb_male_1000_R[ , "L/R" ] <- 0
  
  ukb_male_1000_LR[ , "clean_hv_right" ] <- 0
  ukb_male_1000_R[ , "clean_hv_left" ] <- 0
  
  names(ukb_male_1000_LR)[2] <- "clean_hv"
  names(ukb_male_1000_R)[3] <- "clean_hv"
  
  ukb_male_1000_LR <- ukb_male_1000_LR[ , c( "AGE_Latest" , "clean_hv" , "L/R")]
  ukb_male_1000_R <- ukb_male_1000_R[ , c( "AGE_Latest" , "clean_hv" , "L/R")]
  
  ukb_male_1000_LR <- rbind(ukb_male_1000_LR , ukb_male_1000_R)
  
  bins_sep <- GPR_ANALYSIS(ukb_male_1000)
  bins_com <- GPR_ANALYSIS_TEST_LR(ukb_male_1000_LR)
  
  par(mfrow=c(2,2))
  PLOT_NOMOGRAM(bins_sep , hem = 0 , title = "SEP - LEFT" , xlim = c(47,82) , ylim = c(2700 , 6200) )
  PLOT_NOMOGRAM(bins_sep , hem = 1 , title = "SEP - RIGHT", xlim = c(47,82) , ylim = c(2700 , 6200) )
  PLOT_NOMOGRAM(bins_com , hem = 0 , title = "COM - LEFT" , xlim = c(47,82) , ylim = c(2700 , 6200) )
  PLOT_NOMOGRAM(bins_com , hem = 1 , title = "COM - RIGHT", xlim = c(47,82) , ylim = c(2700 , 6200) )
  
}


#
# Explore run times and data effeciency of SW and GPR methods
#
#
EXPLORE_RUNTIMES <- function(){
  durations_table_sliding_window <- data.frame( train_size=double() , duration=double()) 
  for(num_train_samples in seq(35000,40000,100) ){
    start_time <- Sys.time()
    assign ( paste("bins_male_SW_",num_train_samples,sep="") , WINDOW_ANALYSIS( ukb_img_ex_outliers_male[1:num_train_samples,]) )
    end_time <- Sys.time()
    duration <- as.numeric( end_time - start_time , units="secs")
    print(paste("Trained on:",num_train_samples,"samples. Time taken is:", duration, "Seconds"))
    durations_table_sliding_window[nrow(durations_table_sliding_window) + 1 , ] <- c( num_train_samples , duration )
  }
  
  for(num_train_samples in seq(100,5500,100) )
    assign ( paste("bins_male_gpr","_",num_train_samples,sep="") , read.table(paste( "~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/","bins_male_gpr","_",num_train_samples,sep="" )) )
  durations_table <- read.table("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/durations" )
  
  par(mfrow=c(2,2))
  for( N in c(100,1000,3000,5500)){
    bins = get( paste( "bins_male_gpr_" , N , sep="") )
    PLOT_NOMOGRAM( bins_male , xlim=c(35,90) , ylim=c(3000,6000)  , title = paste("MALE LEFT \n Sliding Window 40k vs GPR trained with",N) )
    lines(bins$mid_age , bins$l.q2.5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q10 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q25 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q50 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q75 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q90 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q95 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q97.5 , lty = 2, lwd = 2 , col = 4)
  }
  par(mfrow=c(2,2))
  for( N in c(100,1000,2000,4000)){
    bins = get( paste( "bins_male_SW_" , N , sep="") )
    PLOT_NOMOGRAM( bins_male , xlim=c(35,90) , ylim=c(3000,6000)  , title = paste("MALE LEFT \n Sliding Window 40k vs Sliding Window trained with",N) )
    lines(bins$mid_age , bins$l.q2.5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q10 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q25 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q50 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q75 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q90 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q95 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q97.5 , lty = 2, lwd = 2 , col = 4)
  }
  
  par(mfrow=c(2,4))
  for( N in c(100,1000,2000,4000)){
    bins = get( paste( "bins_male_SW_" , N , sep="") )
    PLOT_NOMOGRAM( bins_male , xlim=c(35,90) , ylim=c(3000,6000)  , title = paste("MALE LEFT \n Sliding Window 40k vs Sliding Window trained with",N) )
    lines(bins$mid_age , bins$l.q2.5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q10 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q25 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q50 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q75 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q90 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q95 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q97.5 , lty = 2, lwd = 2 , col = 4)
    
    bins = get( paste( "bins_male_gpr_" , N , sep="") )
    PLOT_NOMOGRAM( bins_male , xlim=c(35,90) , ylim=c(3000,6000)  , title = paste("MALE LEFT \n Sliding Window 40k vs GPR trained with",N) )
    lines(bins$mid_age , bins$l.q2.5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q5 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q10 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q25 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q50 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q75 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q90 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q95 , lty = 2, lwd = 2 , col = 4)
    lines(bins$mid_age , bins$l.q97.5 , lty = 2, lwd = 2 , col = 4)
  }
  
}




#
# Explore if using a better GWAS would give more accurate results
#
EXPLORE_GWAS <- function( ukb_male , ukb_female ){
  
  ukb_full <- rbind(ukb_male , ukb_female )
  
  sep <- SEPERATE( ukb_male , 'PRS_TH_1' , 0.3)
  ukb_male_low <- sep$Table_low
  ukb_male_high <- sep$Table_high
  
  sep <- SEPERATE( ukb_female , 'PRS_TH_1' , 0.3)
  ukb_female_low <- sep$Table_low
  ukb_female_high <- sep$Table_high
  
  sep <- SEPERATE( ukb_full , 'PRS_TH_1' , 0.3)
  ukb_full_low <- sep$Table_low
  ukb_full_high <- sep$Table_high
  
  
  bins_male <- WINDOW_ANALYSIS( ukb_male , hv_left_column = "clean_hv_mean")
  bins_female <- WINDOW_ANALYSIS( ukb_female , hv_left_column = "clean_hv_mean")
  
  bins_male_low  <- WINDOW_ANALYSIS( ukb_male_low , hv_left_column = "clean_hv_mean")
  bins_male_high <- WINDOW_ANALYSIS( ukb_male_high , hv_left_column = "clean_hv_mean") 
  
  bins_female_low  <- WINDOW_ANALYSIS( ukb_female_low , hv_left_column = "clean_hv_mean")
  bins_female_high <- WINDOW_ANALYSIS( ukb_female_high , hv_left_column = "clean_hv_mean") 
  
  
  bins_full <- WINDOW_ANALYSIS( ukb_full , hv_left_column = "clean_hv_mean")
  bins_full_high <- WINDOW_ANALYSIS( ukb_full_high , hv_left_column = "clean_hv_mean")
  bins_full_low <- WINDOW_ANALYSIS( ukb_full_low , hv_left_column = "clean_hv_mean")
  
  par(mfrow=c(1,3))
  PLOT_NOMOGRAM_COMPARE(bins_full , bins_full_high )
  PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_high )
  PLOT_NOMOGRAM_COMPARE(bins_female , bins_female_high )
  
  
  summary(lm( clean_hv_mean  ~  PRS_TH_1 + AGE_Latest + (AGE_Latest^2) + SEX + ICV , data = ukb_full))
  summary(lm( clean_hv_mean  ~  PRS_TH_1 + AGE_Latest + (AGE_Latest^2) + SEX + ICV , data = ukb_male))
  summary(lm( clean_hv_mean  ~  PRS_TH_1 + AGE_Latest + (AGE_Latest^2) + SEX + ICV , data = ukb_female))
  
}
