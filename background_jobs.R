library(TTR)
library(smoother)
library(ADNIMERGE)
library(laGP)

source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Nomogram Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/PRS Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Preprocessing Functions.R")

durations_table <- data.frame( traning_size=double() ,  duration_sw=double() ,  diff_sw=double() , duration_gpr=double() , diff_gpr=double() )
durations_table <- read.table("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/durations_compare_2")
bins_male_sw <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male ) 

for(num_train_samples in seq(100,6000,100) ){
  
  start_time <- Sys.time()
  assign ( paste("bins_male_sw","_",num_train_samples,sep="") , WINDOW_ANALYSIS( ukb_img_ex_outliers_male[1:num_train_samples,]) )
  end_time <- Sys.time()
  duration_sw <- as.numeric( end_time - start_time , units="secs")
  
  XX <-matrix( bins_male_sw$mid_age , ncol = 1)
  
  start_time <- Sys.time()
  #assign ( paste("bins_male_gpr_",num_train_samples,sep="") ,
  #         read.table(paste("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/all_compared_parallel",num_train_samples,sep="") ) )
   assign ( paste("bins_male_gpr_",num_train_samples,sep="") , GPR_ANALYSIS( ukb_img_ex_outliers_male[1:num_train_samples,] , XX) )
  end_time <- Sys.time()
  duration_gpr <- as.numeric( end_time - start_time , units="secs")
  
  dif_gpr <- NOMOGRAM_DIFF_INTERPOLATE(  bins_male_sw ,  get(paste("bins_male_gpr","_",num_train_samples,sep="")) )
  dif_sw <- NOMOGRAM_DIFF_INTERPOLATE(  bins_male_sw ,  get(paste("bins_male_sw","_",num_train_samples,sep="")) )
  
  print(paste(num_train_samples,"samples: SW Time:",format(duration_sw, digits = 5),"Seconds. Difference:",format(dif_sw,digits=10),
                                    "mm3. GPR Time:", format(duration_gpr,digits=5) , "Seconds. Difference:",format(dif_gpr,digits=10),"mm3"))
  
  durations_table[nrow(durations_table) + 1 , ] <- c( num_train_samples , duration_sw , dif_sw , duration_gpr  , dif_gpr )
  #temp_durations[nrow(temp_durations) + 1 , ] <- c( num_train_samples , dif_sw )
  
  #write.table( get(paste("bins_male_gpr","_",num_train_samples,sep="")) ,
  #             paste( "~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/","bins_male_gpr_compare_",num_train_samples,sep="" ) )
  #write.table( durations_table , paste( "~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/","durations_compare_2",sep="" ) )
  
  #png(file=paste("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/all_compared_parallel",num_train_samples,".png",sep = ""), width=1500, height=875)
  #par(mfrow=c(2,2))
  #PLOT_NOMOGRAM( get(paste("bins_male_gpr","_",num_train_samples,sep="") ) , title = "GPR Nomogram", xlim = c(53,75) , ylim = c(3100,5100))
  #points( ukb_img_ex_outliers_male[1:num_train_samples,]$AGE_Latest , ukb_img_ex_outliers_male[1:num_train_samples,]$clean_hv_left , col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3))
  #PLOT_NOMOGRAM( get(paste("bins_male_sw","_",num_train_samples,sep="") )  , title = "SW Nomogram", xlim = c(53,75) , ylim = c(3100,5100))
  #points( ukb_img_ex_outliers_male[1:num_train_samples,]$AGE_Latest , ukb_img_ex_outliers_male[1:num_train_samples,]$clean_hv_left , col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3))
  #PLOT_NOMOGRAM_COMPARE( get(paste("bins_male_gpr","_",num_train_samples,sep="") ) ,
  #                       bins_male_sw , left = 1 , col_1 = "red" , col_2 = "blue" ,
  #                       title = paste("GPR Nomogram vs Gold Standard : diff =",floor(dif_gpr),"mm3") ,
  #                       xlim = c(53,75) , ylim = c(3100,5100) )
  #PLOT_NOMOGRAM_COMPARE( get(paste("bins_male_sw","_",num_train_samples,sep="") ) ,
  #                       bins_male_sw , left = 1 , col_1 = "red" , col_2 = "blue" ,
  #                       title = paste("GPR Nomogram vs Gold Standard : diff =",floor(dif_sw),"mm3") ,
  #                       xlim = c(53,75) , ylim = c(3100,5100) )
  #mtext( paste("MALE LEFT NOMOGRAM METHODS WITH",num_train_samples,"TRANING SAMPLES") , side = 3, line = -1.5, outer = TRUE)
  #dev.off()
  }
