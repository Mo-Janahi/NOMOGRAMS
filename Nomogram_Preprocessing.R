# 
# By: Mohammed Janahi


###### MAIN CODE  #######
library(TTR)
library(smoother)
library(ADNIMERGE)
library(laGP)


source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Nomogram Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/PRS Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Preprocessing Functions.R")

# UKB Table is saved as CSV file on local machine
ukb <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200830.csv")

# call the function to filter table, take care of outliers, and correct hippocampal volume for ICV and scan date
x <- PREPROCESS_UKBB(ukb)

# grab the separate gender tables
ukb_img_ex_outliers_male <- x$ukb_img_ex_outliers_male
ukb_img_ex_outliers_female <- x$ukb_img_ex_outliers_female

######### WINDOW ANALYSIS #########

# Now we will generate a table with binned information to begin constructing the nomogram.
bins_male <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male )
bins_female <- WINDOW_ANALYSIS( ukb_img_ex_outliers_female )

XX <- matrix( bins_male$mid_age , ncol=1)
bins_male_gpr <- GPR_ANALYSIS( ukb_img_ex_outliers_male[1:100,] , XX)
######### PLOTTING  NOMOGRAMS ######### 

par(mfrow=c(2,2))
PLOT_NOMOGRAM( bins_male , 1 , "Nomogram:\n male Left HV")
PLOT_NOMOGRAM( bins_male , 0 , "Nomogram:\n male Right HV")
PLOT_NOMOGRAM( bins_female , 1 , "Nomogram:\n Female Left HV")
PLOT_NOMOGRAM( bins_female , 0 , "Nomogram:\n Female Right HV")


###### PRS ANALYSIS #######

MAKE_PRS_PLOTS( ukb_img_ex_outliers_male , ukb_img_ex_outliers_female , c("AGE_Latest","Sex","Genetic.PC.1","Genetic.PC.2","Genetic.PC.3") )


# now we want to select top/bottom 30% of samples and make seperate nomograms and plot them vs old nomograms 

sep <- SEPERATE( ukb_img_ex_outliers_male , 'PRS_TH_1' , 0.3)
ukb_img_ex_outliers_male_low <- sep$Table_low
ukb_img_ex_outliers_male_high <- sep$Table_high

sep <- SEPERATE( ukb_img_ex_outliers_female , 'PRS_TH_1' , 0.3)
ukb_img_ex_outliers_female_low <- sep$Table_low
ukb_img_ex_outliers_female_high <- sep$Table_high

bins_male_low  <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male_low )
bins_male_high <- WINDOW_ANALYSIS( ukb_img_ex_outliers_male_high ) 

bins_female_low  <- WINDOW_ANALYSIS( ukb_img_ex_outliers_female_low )
bins_female_high <- WINDOW_ANALYSIS( ukb_img_ex_outliers_female_high ) 

num_train_samples <- 1000
bins_male_gpr      <- GPR_ANALYSIS( ukb_img_ex_outliers_male[1:num_train_samples,])
bins_male_gpr_low  <- GPR_ANALYSIS( ukb_img_ex_outliers_male_low[1:num_train_samples,] )
bins_male_gpr_high <- GPR_ANALYSIS( ukb_img_ex_outliers_male_high[1:num_train_samples,] ) 

bins_female_gpr      <- GPR_ANALYSIS( ukb_img_ex_outliers_female[1:num_train_samples,] )
bins_female_gpr_low  <- GPR_ANALYSIS( ukb_img_ex_outliers_female_low[1:num_train_samples,] )
bins_female_gpr_high <- GPR_ANALYSIS( ukb_img_ex_outliers_female_high[1:num_train_samples,] ) 

# one plot to show all separate nomograms 
par(mfrow=c(3,4))
for( strat in c("" , "_low" , "_high") )
  for( sex in c("male" , "female"))
    for (hem in c("Left" , "Right")){
      bins <- get(paste("bins_",sex,strat,sep=""))
      title <- paste("HV Nomogram:\n",toupper(sex),hem,substring(strat,2,nchar(strat)),"prs")
      PLOT_NOMOGRAM( bins , (hem=="Left") , title )
    }

# one more plot with nomograms overlaid on each other 
par(mfrow=c(2,4))
for( strat in c("_low" , "_high") )
  for( sex in c("male" , "female"))
    for (hem in c("Left" , "Right")){
      bins_1 <- get(paste("bins_",sex,sep=""))
      bins_2 <- get(paste("bins_",sex,strat,sep=""))
      title <- paste("HV Nomogram:\n",toupper(sex),hem,"Normal vs",substring(strat,2,nchar(strat)),"prs")
      color = "blue" 
      if(strat=="_high") color = "green"
      PLOT_NOMOGRAM_COMPARE( bins_1 , bins_2 , (hem=="Left") , "red" , color , title)
    }


# one more plot with nomograms overlaid on each other 
par(mfrow=c(2,2))
PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL( bins_male , bins_male_low, bins_male_high , 1 , "red" , "blue" , "green", "Male Left \n Normal vs Low vs High prs" )
PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL( bins_male , bins_male_low, bins_male_high , 0 , "red" , "blue" , "green", "Male Right \n Normal vs Low vs High prs" )
PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL( bins_female , bins_female_low , bins_female_high , 1 , "red" , "blue" , "green", "Female Left \n Normal vs Low vs High prs" )
PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL( bins_female , bins_female_low , bins_female_high , 0 , "red" , "blue" , "green", "Female Right \n Normal vs Low vs High prs" )


# this is exploration of where is the best threshold/ percentile to cut samples at 
# THRESHOLD_EXPLORE( ukb_img_ex_outliers_male , ukb_img_ex_outliers_female )


###### NOMOGRAM EVALUATION #########

# we will bring in some samples from ADNI dataset to see if the new nomograms are better than the old ones.

ADNI_table <- merge( adnimerge , ADNIMERGE::ucsffsx51 , by = c("IMAGEUID") )
#ADNI_table <- merge( adnimerge , ADNIMERGE::ucsffsx , by = c("IMAGEUID") ) 
#ADNI_table <- merge( adnimerge , ADNIMERGE::ucsffsx6 , by = c("RID") )

x <- PREPROCESS_ADNI(ADNI_table)

ADNI_filter_male <- x$ADNI_filter_male
ADNI_filter_female <- x$ADNI_filter_female

ADNI_filter_male_AD <- ADNI_filter_male[ (ADNI_filter_male$VISCODE.y=="scmri") & (ADNI_filter_male$DX=="Dementia") , ]
ADNI_filter_male_MCI <- ADNI_filter_male[ (ADNI_filter_male$VISCODE.y=="scmri") & (ADNI_filter_male$DX=="MCI") , ]
ADNI_filter_male_CN <- ADNI_filter_male[ (ADNI_filter_male$VISCODE.y=="scmri") & (ADNI_filter_male$DX=="CN") , ]
ADNI_filter_male_NA <- ADNI_filter_male[ (ADNI_filter_male$VISCODE.y=="scmri") & (is.na(ADNI_filter_male$DX)) , ]


ADNI_filter_female_AD <- ADNI_filter_female[ (ADNI_filter_female$VISCODE.y=="scmri") & (ADNI_filter_female$DX=="Dementia") , ]
ADNI_filter_female_MCI <- ADNI_filter_female[ (ADNI_filter_female$VISCODE.y=="scmri") & (ADNI_filter_female$DX=="MCI") , ]
ADNI_filter_female_CN <- ADNI_filter_female[ (ADNI_filter_female$VISCODE.y=="scmri") & (ADNI_filter_female$DX=="CN") , ]
ADNI_filter_female_NA <- ADNI_filter_female[ (ADNI_filter_female$VISCODE.y=="scmri") & (is.na(ADNI_filter_female$DX)) , ]


par(mfrow=c(2,2))
PLOT_NOMOGRAM( bins_male , 1 , paste("Nomogram:\n Male Right HV") , xlim = c(52,76))
points((ADNI_filter_male_AD$TRUE.AGE) , (ADNI_filter_male_AD$clean_hv_right) , pch=17 , col="red")
points((ADNI_filter_male_MCI$TRUE.AGE) , (ADNI_filter_male_MCI$clean_hv_right) , pch=16 , col="green")
points((ADNI_filter_male_CN$TRUE.AGE) , (ADNI_filter_male_CN$clean_hv_right) , pch=15 , col="blue")
points((ADNI_filter_male_NA$TRUE.AGE) , (ADNI_filter_male_NA$clean_hv_right) , pch=15 , col="grey")

PLOT_NOMOGRAM( bins_male , 0 , paste("Nomogram:\n Male Left HV") , xlim = c(52,76))
points((ADNI_filter_male_AD$TRUE.AGE) , (ADNI_filter_male_AD$clean_hv_left) , pch=17 , col="red")
points((ADNI_filter_male_MCI$TRUE.AGE) , (ADNI_filter_male_MCI$clean_hv_left) , pch=16 , col="green")
points((ADNI_filter_male_CN$TRUE.AGE) , (ADNI_filter_male_CN$clean_hv_left) , pch=15 , col="blue")
points((ADNI_filter_male_NA$TRUE.AGE) , (ADNI_filter_male_NA$clean_hv_left) , pch=15 , col="grey")

PLOT_NOMOGRAM( bins_female , 1 , paste("Nomogram:\n Female Right HV") , xlim = c(52,76))
points((ADNI_filter_female_AD$TRUE.AGE) , (ADNI_filter_female_AD$clean_hv_right) , pch=17 , col="red")
points((ADNI_filter_female_MCI$TRUE.AGE) , (ADNI_filter_female_MCI$clean_hv_right) , pch=16 , col="green")
points((ADNI_filter_female_CN$TRUE.AGE) , (ADNI_filter_female_CN$clean_hv_right) , pch=15 , col="blue")
points((ADNI_filter_female_NA$TRUE.AGE) , (ADNI_filter_female_NA$clean_hv_right) , pch=15 , col="grey")

PLOT_NOMOGRAM( bins_female , 0 , paste("Nomogram:\n Female Left HV") , xlim = c(52,76))
points((ADNI_filter_female_AD$TRUE.AGE) , (ADNI_filter_female_AD$clean_hv_left) , pch=17 , col="red")
points((ADNI_filter_female_MCI$TRUE.AGE) , (ADNI_filter_female_MCI$clean_hv_left) , pch=16 , col="green")
points((ADNI_filter_female_CN$TRUE.AGE) , (ADNI_filter_female_CN$clean_hv_left) , pch=15 , col="blue")
points((ADNI_filter_female_NA$TRUE.AGE) , (ADNI_filter_female_NA$clean_hv_left) , pch=15 , col="grey")



MAKE_PRS_PLOTS( ADNI_filter_male , ADNI_filter_female , c("TRUE.AGE","Sex","PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10") ,  
                hv_l_column = "clean_hv_left", hv_r_column = "clean_hv_right",
                hv_bl_column = "clean_hv_bilateral" )


# create a bunch of variables with a loop rather than write them all out
# variables are in the form of:
# M_L_AD <- AVERAGE_PERCENTILE( ADNI_filter_male_AD , bins_male , "clean_hv_left" , 1 )
# and we are generating:  M_L_AD, M_L_MCI, M_L_CN, M_R_AD, M_R_MCI, M_R_CN, F_L_AD, F_L_MCI, F_L_CN, F_R_AD, F_R_MCI, F_R_CN
for( x in c("male","female") )
  for( y in c("left","right") )
    for(z in c("AD","MCI","CN","NA") )
      assign( paste(toupper(substring(x,1,1)),"_",toupper(substring(y,1,1)),"_",z,sep = "") ,
              AVERAGE_PERCENTILE( get(paste("ADNI_filter_",x,"_",z,sep = "")) , 
                                  get(paste("bins_",x,sep = "")) , 
                                  paste("clean_hv_",y,sep = "") , 
                                  (y=="left") , xlim = NA ) )



par(mfrow=c(1,1))
boxplot(  c(M_L_AD , M_L_MCI , M_L_CN , M_R_AD , M_R_MCI , M_R_CN , F_L_AD , F_L_MCI , F_L_CN , F_R_AD , F_R_MCI , F_R_CN ) ~ 
            c( (M_L_AD * 0 + 1),(M_L_MCI * 0 + 2), (M_L_CN * 0 + 3),(M_R_AD * 0 + 4), (M_R_MCI * 0 + 5),(M_R_CN * 0 + 6),
               (F_L_AD * 0 + 7),(F_L_MCI * 0 + 8), (F_L_CN * 0 + 9),(F_R_AD * 0 + 10), (F_R_MCI * 0 + 11),(F_R_CN * 0 + 12)  ), 
          col=c( rep(rgb(0,0,1,0.2) , 3), rep(rgb(0,0,1,0.6) , 3), 
                 rep(rgb(1,0,0,0.6) , 3), rep(rgb(1,0,0,0.9) , 3) ),
          outline=FALSE , ylab = "Percentile" , xlab = "Strata of Data" , main = "HV Percentile distribution in GPR-Method Nomogram")

abline(h=50 , col='grey' , lty=1)



sep <- SEPERATE( ADNI_filter_male , 'PRS_TH_1' , 0.5)
ADNI_filter_male_low <- sep$Table_low
ADNI_filter_male_high <- sep$Table_high

sep <- SEPERATE( ADNI_filter_female , 'PRS_TH_1' , 0.5)
ADNI_filter_female_low <- sep$Table_low
ADNI_filter_female_high <- sep$Table_high

ADNI_filter_male_high_AD <- ADNI_filter_male_high[ (ADNI_filter_male_high$VISCODE.y=="scmri") & (ADNI_filter_male_high$DX=="Dementia") , ]
ADNI_filter_male_high_MCI <- ADNI_filter_male_high[ (ADNI_filter_male_high$VISCODE.y=="scmri") & (ADNI_filter_male_high$DX=="MCI") , ]
ADNI_filter_male_high_CN <- ADNI_filter_male_high[ (ADNI_filter_male_high$VISCODE.y=="scmri") & (ADNI_filter_male_high$DX=="CN") , ]
ADNI_filter_male_high_NA <- ADNI_filter_male_high[ (ADNI_filter_male_high$VISCODE.y=="scmri") & (is.na(ADNI_filter_male_high$DX)) , ]

ADNI_filter_male_low_AD <- ADNI_filter_male_low[ (ADNI_filter_male_low$VISCODE.y=="scmri") & (ADNI_filter_male_low$DX=="Dementia") , ]
ADNI_filter_male_low_MCI <- ADNI_filter_male_low[ (ADNI_filter_male_low$VISCODE.y=="scmri") & (ADNI_filter_male_low$DX=="MCI") , ]
ADNI_filter_male_low_CN <- ADNI_filter_male_low[ (ADNI_filter_male_low$VISCODE.y=="scmri") & (ADNI_filter_male_low$DX=="CN") , ]
ADNI_filter_male_low_NA <- ADNI_filter_male_low[ (ADNI_filter_male_low$VISCODE.y=="scmri") & (is.na(ADNI_filter_male_low$DX)) , ]

ADNI_filter_female_high_AD <- ADNI_filter_female_high[ (ADNI_filter_female_high$VISCODE.y=="scmri") & (ADNI_filter_female_high$DX=="Dementia") , ]
ADNI_filter_female_high_MCI <- ADNI_filter_female_high[ (ADNI_filter_female_high$VISCODE.y=="scmri") & (ADNI_filter_female_high$DX=="MCI") , ]
ADNI_filter_female_high_CN <- ADNI_filter_female_high[ (ADNI_filter_female_high$VISCODE.y=="scmri") & (ADNI_filter_female_high$DX=="CN") , ]
ADNI_filter_female_high_NA <- ADNI_filter_female_high[ (ADNI_filter_female_high$VISCODE.y=="scmri") & (is.na(ADNI_filter_female_high$DX)) , ]

ADNI_filter_female_low_AD <- ADNI_filter_female_low[ (ADNI_filter_female_low$VISCODE.y=="scmri") & (ADNI_filter_female_low$DX=="Dementia") , ]
ADNI_filter_female_low_MCI <- ADNI_filter_female_low[ (ADNI_filter_female_low$VISCODE.y=="scmri") & (ADNI_filter_female_low$DX=="MCI") , ]
ADNI_filter_female_low_CN <- ADNI_filter_female_low[ (ADNI_filter_female_low$VISCODE.y=="scmri") & (ADNI_filter_female_low$DX=="CN") , ]
ADNI_filter_female_low_NA <- ADNI_filter_female_low[ (ADNI_filter_female_low$VISCODE.y=="scmri") & (is.na(ADNI_filter_female_low$DX)) , ]

# create a bunch of variables with a loop rather than write them all out
# variables are in the form of:
#M_L_CN_HN <- AVERAGE_PERCENTILE( ADNI_filter_male_high_CN , bins_male , "clean_hv_left" , 1 )
# and we are generating:
#  M_L_CN_HN, M_L_CN_H, M_L_CN_LN, M_L_CN_L, M_R_CN_HN, M_R_CN_H, M_R_CN_LN, M_R_CN_L,
#  F_L_CN_HN, F_L_CN_H, F_L_CN_LN, F_L_CN_L, F_R_CN_HN, F_R_CN_H, F_R_CN_LN, F_R_CN_L 
for( x in c("male","female") )
  for( y in c("left","right") )
    for(z in c("AD","MCI","CN","NA") )
      for(j in c("high","low")){
        assign( paste(toupper(substring(x,1,1)),"_",toupper(substring(y,1,1)),"_",z,"_",toupper(substring(j,1,1)),"N",sep = "") ,
                AVERAGE_PERCENTILE( get(paste("ADNI_filter_",x,"_",j,"_",z,sep = "")) , 
                                    get(paste("bins_",x,"_gpr_accum",sep = "")) , 
                                    paste("clean_hv_",y,sep = "") , (y=="left") ) )
        assign( paste(toupper(substring(x,1,1)),"_",toupper(substring(y,1,1)),"_",z,"_",toupper(substring(j,1,1)),sep = "") ,
                AVERAGE_PERCENTILE( get(paste("ADNI_filter_",x,"_",j,"_",z,sep = "")) , 
                                    get(paste("bins_",x,"_gpr_accum_",j,sep = "")) , 
                                    paste("clean_hv_",y,sep = "") , (y=="left") ) )
      }



par(mfrow=c(1,1))

boxplot(  c(M_L_CN_HN , M_L_CN_H , M_L_CN_LN , M_L_CN_L , M_R_CN_HN , M_R_CN_H , M_R_CN_LN , M_R_CN_L ,
            F_L_CN_HN , F_L_CN_H , F_L_CN_LN , F_L_CN_L , F_R_CN_HN , F_R_CN_H , F_R_CN_LN , F_R_CN_L ) ~ 
            c( (M_L_CN_HN*0 + 1) , (M_L_CN_H*0 + 2), (M_L_CN_LN*0 + 3),(M_L_CN_L*0 + 4),
               (M_R_CN_HN*0 + 5) , (M_R_CN_H*0 + 6),(M_R_CN_LN*0 + 7),(M_R_CN_L*0 + 8),
               (F_L_CN_HN*0 + 9) , (F_L_CN_H*0 + 10),(F_L_CN_LN*0 + 11),(F_L_CN_L*0 + 12),
               (F_R_CN_HN*0 + 13), (F_R_CN_H*0 + 14),(F_R_CN_LN*0 + 15),(F_R_CN_L*0 + 16) ),
          col= c( rep(c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2) ), 4) ),
          outline=FALSE, ylab = "Percentile" , xlab = "Strata of Data" , main = "CN HV Percentile distribution in GPR-Method Nomogram + ADNI CN"
)
abline(h=50 , col='grey' , lty=1)

boxplot(  c(M_L_MCI_HN , M_L_MCI_H , M_L_MCI_LN , M_L_MCI_L ,
            M_R_MCI_HN , M_R_MCI_H , M_R_MCI_LN , M_R_MCI_L ,
            F_L_MCI_HN , F_L_MCI_H , F_L_MCI_LN , F_L_MCI_L ,
            F_R_MCI_HN , F_R_MCI_H , F_R_MCI_LN , F_R_MCI_L ) ~ 
            c( (M_L_MCI_HN*0 + 1),(M_L_MCI_H*0 + 2),(M_L_MCI_LN*0 + 3),(M_L_MCI_L*0 + 4),
               (M_R_MCI_HN*0 + 5),(M_R_MCI_H*0 + 6),(M_R_MCI_LN*0 + 7),(M_R_MCI_L*0 + 8),
               (F_L_MCI_HN*0 + 9),(F_L_MCI_H*0 + 10),(F_L_MCI_LN*0 + 11),(F_L_MCI_L*0 + 12),
               (F_R_MCI_HN*0 + 13),(F_R_MCI_H*0 + 14),(F_R_MCI_LN*0 + 15),(F_R_MCI_L*0 + 16) ) ,
          col= c( rep(c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2) ), 4) ),
          outline=FALSE , ylab = "Percentile" , xlab = "Strata of Data" 
)

boxplot(  c(M_L_AD_HN , M_L_AD_H , M_L_AD_LN , M_L_AD_L ,
            M_R_AD_HN , M_R_AD_H , M_R_AD_LN , M_R_AD_L ,
            F_L_AD_HN , F_L_AD_H , F_L_AD_LN , F_L_AD_L ,
            F_R_AD_HN , F_R_AD_H , F_R_AD_LN , F_R_AD_L ) ~ 
            c( (1:length(M_L_AD_HN) * 0 + 1),(1:length(M_L_AD_H) * 0 + 2),(1:length(M_L_AD_LN) * 0 + 3),(1:length(M_L_AD_L) * 0 + 4),
               (1:length(M_R_AD_HN) * 0 + 5),(1:length(M_R_AD_H) * 0 + 6),(1:length(M_R_AD_LN) * 0 + 7),(1:length(M_R_AD_L) * 0 + 8),
               (1:length(F_L_AD_HN) * 0 + 9),(1:length(F_L_AD_H) * 0 + 10),(1:length(F_L_AD_LN) * 0 + 11),(1:length(F_L_AD_L) * 0 + 12),
               (1:length(F_R_AD_HN) * 0 + 13),(1:length(F_R_AD_H) * 0 + 14),(1:length(F_R_AD_LN) * 0 + 15),(1:length(F_R_AD_L) * 0 + 16) ) ,
          col= c( rep(c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2) ), 4) ),
          outline=FALSE , ylab = "Percentile" , xlab = "Strata of Data" , main = "AD HV Percentile distribution in GPR-Method Nomogram + ADNI CN"
)

par(mfrow=c(3,1))

PLOT_NOMOGRAM( bins_female , 0 , paste("Nomogram:\n Female Left HV"))
points((ADNI_filter_female_AD$TRUE.AGE) , (ADNI_filter_female_AD$clean_hv_left) , pch=17 , col="red")
points((ADNI_filter_female_MCI$TRUE.AGE) , (ADNI_filter_female_MCI$clean_hv_left) , pch=16 , col="green")
points((ADNI_filter_female_CN$TRUE.AGE) , (ADNI_filter_female_CN$clean_hv_left) , pch=15 , col="blue")
points((ADNI_filter_female_NA$TRUE.AGE) , (ADNI_filter_female_NA$clean_hv_left) , pch=15 , col="grey")

PLOT_NOMOGRAM( bins_female_high , 0 , paste("Nomogram:\n Female Left HV"))
points((ADNI_filter_female_high_AD$TRUE.AGE) , (ADNI_filter_female_high_AD$clean_hv_left) , pch=17 , col="red")
points((ADNI_filter_female_high_MCI$TRUE.AGE) , (ADNI_filter_female_high_MCI$clean_hv_left) , pch=16 , col="green")
points((ADNI_filter_female_high_CN$TRUE.AGE) , (ADNI_filter_female_high_CN$clean_hv_left) , pch=15 , col="blue")
points((ADNI_filter_female_high_NA$TRUE.AGE) , (ADNI_filter_female_high_NA$clean_hv_left) , pch=15 , col="grey")

PLOT_NOMOGRAM( bins_female_low , 0 , paste("Nomogram:\n Female Left HV"))
points((ADNI_filter_female_low_AD$TRUE.AGE) , (ADNI_filter_female_low_AD$clean_hv_left) , pch=17 , col="red")
points((ADNI_filter_female_low_MCI$TRUE.AGE) , (ADNI_filter_female_low_MCI$clean_hv_left) , pch=16 , col="green")
points((ADNI_filter_female_low_CN$TRUE.AGE) , (ADNI_filter_female_low_CN$clean_hv_left) , pch=15 , col="blue")
points((ADNI_filter_female_low_NA$TRUE.AGE) , (ADNI_filter_female_low_NA$clean_hv_left) , pch=15 , col="grey")


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


data <- rbind( ukb_img_ex_outliers_male[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
               setNames( ADNI_filter_male_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                         c("AGE_Latest","clean_hv_left","clean_hv_right")) )
bins_male_gpr_accum <- GPR_ANALYSIS(data)

data <- rbind( ukb_img_ex_outliers_female[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
               setNames( ADNI_filter_male_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                         c("AGE_Latest","clean_hv_left","clean_hv_right")) )
bins_female_gpr_accum <- GPR_ANALYSIS(data)

data <- rbind( ukb_img_ex_outliers_male_high[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
               setNames( ADNI_filter_male_high_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                         c("AGE_Latest","clean_hv_left","clean_hv_right")) )
bins_male_gpr_accum_high <- GPR_ANALYSIS(data)

data <- rbind( ukb_img_ex_outliers_female_high[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
               setNames( ADNI_filter_male_high_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                         c("AGE_Latest","clean_hv_left","clean_hv_right")) )
bins_female_gpr_accum_high <- GPR_ANALYSIS(data)

data <- rbind( ukb_img_ex_outliers_male_low[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
               setNames( ADNI_filter_male_low_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                         c("AGE_Latest","clean_hv_left","clean_hv_right")) )
bins_male_gpr_accum_low <- GPR_ANALYSIS(data)

data <- rbind( ukb_img_ex_outliers_female_low[1:1000 , c("AGE_Latest","clean_hv_left","clean_hv_right")] , 
               setNames( ADNI_filter_male_low_CN[ , c("TRUE.AGE","clean_hv_left","clean_hv_right")] , 
                         c("AGE_Latest","clean_hv_left","clean_hv_right")) )
bins_female_gpr_accum_low <- GPR_ANALYSIS(data)


par(mfrow=c(1,2))
PLOT_NOMOGRAM( bins_4 , 1 , paste("Nomogram:\n Male Right HV"))
points((ADNI_filter_male_AD$TRUE.AGE) , (ADNI_filter_male_AD$clean_hv_right) , pch=17 , col="red")
points((ADNI_filter_male_MCI$TRUE.AGE) , (ADNI_filter_male_MCI$clean_hv_right) , pch=16 , col="green")
points((ADNI_filter_male_CN$TRUE.AGE) , (ADNI_filter_male_CN$clean_hv_right) , pch=15 , col="blue")
points((ADNI_filter_male_NA$TRUE.AGE) , (ADNI_filter_male_NA$clean_hv_right) , pch=15 , col="grey")

PLOT_NOMOGRAM( bins_4 , 0 , paste("Nomogram:\n Male Left HV"))
points((ADNI_filter_male_AD$TRUE.AGE) , (ADNI_filter_male_AD$clean_hv_left) , pch=17 , col="red")
points((ADNI_filter_male_MCI$TRUE.AGE) , (ADNI_filter_male_MCI$clean_hv_left) , pch=16 , col="green")
points((ADNI_filter_male_CN$TRUE.AGE) , (ADNI_filter_male_CN$clean_hv_left) , pch=15 , col="blue")
points((ADNI_filter_male_NA$TRUE.AGE) , (ADNI_filter_male_NA$clean_hv_left) , pch=15 , col="grey")

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
##### LONGITUDINAL ANALYSIS #########

par(mfrow=c(1,1))
PLOT_NOMOGRAM( bins_female , 0 , paste("Nomogram:\n Female Left HV"))
table <- ADNI_filter_female[ c("RID.x" , "TRUE.AGE" , "clean_hv_left", "clean_hv_right" , "clean_hv_bilateral" , "DX")]
hv_col <- "clean_hv_left"
IDS <- unique(table$RID.x)
for (ID in IDS ){
  p_table <- table[table$RID.x==ID , ]
  p_table <- p_table[ order(p_table$TRUE.AGE) , ]
  if( length(p_table[ , hv_col]) > 1)
    if(p_table[1,"TRUE.AGE"] < 73 )
    if( sum(!is.na(unique(p_table$DX))) > 1 ){
      lines( p_table$TRUE.AGE , p_table[,hv_col] , col=rgb(0,0,0, alpha = 0.5))
      points((p_table[ p_table$DX=="Dementia","TRUE.AGE"]) , (p_table[ p_table$DX=="Dementia",hv_col]) , pch=17 , col="red")
      points((p_table[ p_table$DX=="MCI","TRUE.AGE"]) , (p_table[ p_table$DX=="MCI",hv_col]) , pch=16 , col="green")
      points((p_table[ p_table$DX=="CN","TRUE.AGE"]) , (p_table[ p_table$DX=="CN",hv_col]), pch=15 , col="blue")
      points((p_table[ is.na(p_table$DX),"TRUE.AGE"]) , (p_table[ is.na(p_table$DX),hv_col]) , pch=15 , col="grey")
    }
}
