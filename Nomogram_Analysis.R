# 
# By: Mohammed Janahi


###### MAIN CODE  #######
library(TTR)
library(smoother)
library(ADNIMERGE)
library(laGP)
library(mvtnorm)
library(dplyr)


source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Nomogram Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/PRS Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Preprocessing Functions.R")
source("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/Exploration Functions.R")

# UKB Table is saved as CSV file on local machine
ukb <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200830.csv")

# call the function to filter table, take care of outliers, and correct hippocampal volume for ICV and scan date
x <- PREPROCESS_UKBB(ukb)

# grab the separate gender tables
ukb_male <- x$ukb_img_ex_outliers_male
ukb_female <- x$ukb_img_ex_outliers_female

######### WINDOW ANALYSIS #########

# Now we will generate a table with binned information to begin constructing the nomogram.
bins_male <- WINDOW_ANALYSIS( ukb_male )
bins_female <- WINDOW_ANALYSIS( ukb_female )

# Loading pre-trained GP models
gp_model_male_mean  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_MEAN_HV_MALE")
gp_model_male_left  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_LEFT_HV_MALE")
gp_model_male_right <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_RIGHT_HV_MALE")

# or re-train models 
#gp_model_male_mean <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_bilateral")
#gp_model_male_left <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_left")
#gp_model_male_right <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_right")

ages_male <- seq( min(ukb_male$AGE_Latest) , max(ukb_male$AGE_Latest) , 0.25)
bins_male_gpr <- GPR_MODEL_TO_BINS(gp_model_male_mean , gp_model_male_left , gp_model_male_right , ages_male)

gp_model_female_mean  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_MEAN_HV_FEMALE")
gp_model_female_left  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_LEFT_HV_FEMALE")
gp_model_female_right <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_RIGHT_HV_FEMALE")

#gp_model_female_mean <- GPR_MODEL(ukb_female[1:1000,])
#gp_model_female_left <- GPR_MODEL(ukb_female[1:1000,])
#gp_model_female_right <- GPR_MODEL(ukb_female[1:1000,])

ages_female <- seq( min(ukb_female$AGE_Latest) , max(ukb_female$AGE_Latest) , 0.25)
bins_female_gpr <- GPR_MODEL_TO_BINS(gp_model_female_mean , gp_model_female_left , gp_model_female_right , ages_female)

######### PLOTTING  NOMOGRAMS ######### 

# Figure 1 : 
par(mfrow=c(2,2))
PLOT_NOMOGRAM( bins_male_gpr , "left" , "Male Left" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_male_gpr , "right" , "Male Right", ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female_gpr , "left" , "Female Left", ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female_gpr , "right" , "Female Right", ylim = c(2700,5300) , xlim = c(40,85))


# Figure 2 : 
par(mfrow=c(1,2))
PLOT_NOMOGRAM( bins_male , "bilateral" , "Sliding Window Method" , ylim = c(1500,5300) , xlim = c(40,95))
points(ADNI_male[ADNI_male$TRUE.AGE <= 73, "TRUE.AGE"] , ADNI_male[ADNI_male$TRUE.AGE<=73 , "clean_hv_bilateral"] , col = "grey" , pch=20)
points(ADNI_male[ADNI_male$TRUE.AGE > 73, "TRUE.AGE"] , ADNI_male[ADNI_male$TRUE.AGE>73 , "clean_hv_bilateral"] , col = "red" , pch=20)
PLOT_NOMOGRAM( bins_male_gpr , "bilateral" , "Gaussian Process Regresssion", ylim = c(1500,5300) , xlim = c(40,95))
points(ADNI_male[ADNI_male$TRUE.AGE <= 82, "TRUE.AGE"] , ADNI_male[ADNI_male$TRUE.AGE<=82 , "clean_hv_bilateral"] , col = "grey", pch=20)
points(ADNI_male[ADNI_male$TRUE.AGE > 82, "TRUE.AGE"] , ADNI_male[ADNI_male$TRUE.AGE>82 , "clean_hv_bilateral"] , col = "red", pch=20)

# Supplementary Figure 1 :
par(mfrow=c(3,4))
PLOT_NOMOGRAM( bins_male , "left" , "Male Left SWM" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_male_gpr , "left" , "Male Left GPR" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female , "left" , "Female Left SWM" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female_gpr , "left" , "Female Left GPR", ylim = c(2700,5300) , xlim = c(40,85))

PLOT_NOMOGRAM( bins_male , "right" , "Male Right SWM" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_male_gpr , "right" , "Male Right GPR" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female , "right" , "Female Right SWM" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female_gpr , "right" , "Female Right GPR", ylim = c(2700,5300) , xlim = c(40,85))

PLOT_NOMOGRAM( bins_male , "bilateral" , "Male Bilateral SWM" )
PLOT_NOMOGRAM( bins_male_gpr , "bilateral" , "Male Bilateral GPR" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female , "bilateral" , "Female Bilateral SWM" , ylim = c(2700,5300) , xlim = c(40,85))
PLOT_NOMOGRAM( bins_female_gpr , "bilateral" , "Female Bilateral GPR", ylim = c(2700,5300) , xlim = c(40,85))


###### PRS ANALYSIS #######

MAKE_PRS_PLOTS( ukb_male , ukb_female , 
                confounder_columns = c("AGE_Latest","Sex","Genetic.PC.1","Genetic.PC.2","Genetic.PC.3",
                                       "Genetic.PC.4","Genetic.PC.5","Genetic.PC.6","Genetic.PC.7",
                                       "Genetic.PC.8","Genetic.PC.9","Genetic.PC.10") )


# select top/bottom 30% of samples and make separate nomograms and plot them vs nomograms based on full sets

sep <- SEPERATE( ukb_male , 'PRS_TH_1' , 0.3)
ukb_male_low <- sep$Table_low
ukb_male_high <- sep$Table_high

sep <- SEPERATE( ukb_female , 'PRS_TH_1' , 0.3)
ukb_female_low <- sep$Table_low
ukb_female_high <- sep$Table_high

bins_male_low  <- WINDOW_ANALYSIS( ukb_male_low )
bins_male_high <- WINDOW_ANALYSIS( ukb_male_high ) 

bins_female_low  <- WINDOW_ANALYSIS( ukb_female_low )
bins_female_high <- WINDOW_ANALYSIS( ukb_female_high ) 

#bins_male_gpr_low  <- GPR_ANALYSIS( ukb_male_low[1:1000,] )
#bins_male_gpr_high <- GPR_ANALYSIS( ukb_male_high[1:1000,] ) 

#bins_female_gpr_low  <- GPR_ANALYSIS( ukb_female_low[1:1000,] )
#bins_female_gpr_high <- GPR_ANALYSIS( ukb_female_high[1:1000,] ) 


# Loading pre-trained GP models
gp_model_male_mean_high  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_MEAN_HV_MALE_HIGH")
gp_model_male_left_high  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_LEFT_HV_MALE_HIGH")
gp_model_male_right_high <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_RIGHT_HV_MALE_HIGH")

# or re-train models 
#gp_model_male_mean <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_bilateral")
#gp_model_male_left <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_left")
#gp_model_male_right <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_right")

ages_male <- seq( min(ukb_male$AGE_Latest) , max(ukb_male$AGE_Latest) , 0.25)
bins_male_gpr_high <- GPR_MODEL_TO_BINS(gp_model_male_mean_high , gp_model_male_left_high , gp_model_male_right_high , ages_male)

gp_model_female_mean_high  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_MEAN_HV_FEMALE_HIGH")
gp_model_female_left_high  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_LEFT_HV_FEMALE_HIGH")
gp_model_female_right_high <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_RIGHT_HV_FEMALE_HIGH")

#gp_model_female_mean <- GPR_MODEL(ukb_female[1:1000,])
#gp_model_female_left <- GPR_MODEL(ukb_female[1:1000,])
#gp_model_female_right <- GPR_MODEL(ukb_female[1:1000,])

ages_female <- seq( min(ukb_female$AGE_Latest) , max(ukb_female$AGE_Latest) , 0.25)
bins_female_gpr_high <- GPR_MODEL_TO_BINS(gp_model_female_mean_high , gp_model_female_left_high , gp_model_female_right_high , ages_female)

# Loading pre-trained GP models
gp_model_male_mean_low  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_MEAN_HV_MALE_LOW")
gp_model_male_left_low  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_LEFT_HV_MALE_LOW")
gp_model_male_right_low <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_RIGHT_HV_MALE_LOW")

# or re-train models 
#gp_model_male_mean <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_bilateral")
#gp_model_male_left <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_left")
#gp_model_male_right <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_right")

ages_male <- seq( min(ukb_male$AGE_Latest) , max(ukb_male$AGE_Latest) , 0.25)
bins_male_gpr_low <- GPR_MODEL_TO_BINS(gp_model_male_mean_low , gp_model_male_left_low , gp_model_male_right_low , ages_male)

gp_model_female_mean_low  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_MEAN_HV_FEMALE_LOW")
gp_model_female_left_low  <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_LEFT_HV_FEMALE_LOW")
gp_model_female_right_low <- LOAD_GP("~/Desktop/UKB NOMOGRAM PROJECT/Scripts/NOMOGRAMS/TEMP_GPR_MODELS/FINAL_GP_RIGHT_HV_FEMALE_LOW")

#gp_model_female_mean <- GPR_MODEL(ukb_female[1:1000,])
#gp_model_female_left <- GPR_MODEL(ukb_female[1:1000,])
#gp_model_female_right <- GPR_MODEL(ukb_female[1:1000,])

ages_female <- seq( min(ukb_female$AGE_Latest) , max(ukb_female$AGE_Latest) , 0.25)
bins_female_gpr_low <- GPR_MODEL_TO_BINS(gp_model_female_mean_low , gp_model_female_left_low , gp_model_female_right_low , ages_female)



# plotting is based on the variables above. to avoid writing many lines, 
# I loop through variable names, and get the value of the variable in each loop

# one plot to show all separate nomograms 
par(mfrow=c(2,3))
for( sex in c("male" , "female"))
  for( strat in c("" , "_low" , "_high") ){
      bins <- get(paste("bins_",sex,"_gpr",strat,sep=""))
      title <- paste("HV",toupper(sex),strat,"prs")
      PLOT_NOMOGRAM( bins , hem="bilateral" , title=title )
  }

# one more plot with nomograms overlaid on each other 
par(mfrow=c(3,6))
for( hem in c("left","right","bilateral") )
for( sex in c("male" , "female")){
  bins_1 <- get(paste("bins",sex,"gpr",sep="_"))
  PLOT_NOMOGRAM( bins_1 , title=paste(toupper(sex),toupper(hem),"NOMOGRAM") , hem = hem )
  for( strat in c("high" , "low") ){
      bins_2 <- get(paste("bins",sex,"gpr",strat,sep="_"))
      title <- paste(toupper(strat),"PRS")
      PLOT_NOMOGRAM_COMPARE( bins_1 , bins_2 , hem = hem , "black" , "red" , title )
    }
}
# this last set of plots shows what the plotting code above would look like without the loops

# one more plot with nomograms overlaid on each other 
par(mfrow=c(1,2))
PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL( bins_male_gpr , bins_male_gpr_low, bins_male_gpr_high ,
                                     "bilateral" , "black" , "red" , "blue", "Male Normal vs Low vs High prs" ,
                                     shade = TRUE, xlim=c(46,80))

PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL( bins_female_gpr , bins_female_gpr_low , bins_female_gpr_high ,
                                     "bilateral" , "black" , "red" , "blue", "Female Normal vs Low vs High prs",
                                     shade = TRUE, xlim=c(46,80) )


###### NOMOGRAM EVALUATION #########

# we will bring in some samples from ADNI dataset to see if the new nomograms are better than the old ones.

ADNI_table <- merge( adnimerge , ADNIMERGE::ucsffsx51 , by = c("IMAGEUID") )
#ADNI_table <- merge( adnimerge , ADNIMERGE::ucsffsx , by = c("IMAGEUID") ) 
#ADNI_table <- merge( adnimerge , ADNIMERGE::ucsffsx6 , by = c("RID") )

x <- PREPROCESS_ADNI(ADNI_table)

ADNI_male <- x$ADNI_male
ADNI_male_AD <- filter(ADNI_male , VISCODE.y=="scmri" & DX=="Dementia" )
ADNI_male_MCI <- filter(ADNI_male , VISCODE.y=="scmri" & DX=="MCI" ) 
ADNI_male_CN <- filter(ADNI_male , VISCODE.y=="scmri" & DX=="CN" ) 
ADNI_male_NA <- filter(ADNI_male , VISCODE.y=="scmri" & is.na(DX) ) 

ADNI_female <- x$ADNI_female
ADNI_female_AD <- filter(ADNI_female , VISCODE.y=="scmri" & DX=="Dementia" )
ADNI_female_MCI <- filter(ADNI_female , VISCODE.y=="scmri" & DX=="MCI" )
ADNI_female_CN <- filter(ADNI_female , VISCODE.y=="scmri" & DX=="CN" ) 
ADNI_female_NA <- filter(ADNI_female , VISCODE.y=="scmri" & is.na(DX) )


par(mfrow=c(1,2))
PLOT_NOMOGRAM_ADNI(bins_male_gpr , ADNI_male , hem="bilateral" , title=paste("Nomogram Male HV") , ylim = c(1400,5300) , xlim=c(45,100))
PLOT_NOMOGRAM_ADNI(bins_female_gpr , ADNI_female , hem="bilateral" , title=paste("Nomogram Female HV") , ylim = c(1400,5300), xlim=c(45,100))



MAKE_PRS_PLOTS( ADNI_male , ADNI_female , c("TRUE.AGE","Sex","PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10") ,  
                hv_l_column = "clean_hv_left", hv_r_column = "clean_hv_right",
                hv_bl_column = "clean_hv_bilateral" )


# create a bunch of variables with a loop rather than write them all out
# variables are in the form of:
# M_AD <- AVERAGE_PERCENTILE( ADNI_male_AD , bins_male , "clean_hv_left" , 1 )
# and we are generating:  M_AD, M_MCI, M_CN, M_R_AD, M_R_MCI, M_R_CN, F_AD, F_MCI, F_CN, F_R_AD, F_R_MCI, F_R_CN
for( x in c("male","female") )
    for(z in c("AD","MCI","CN","NA") )
      assign( paste(toupper(substring(x,1,1)),"_",z,sep = "") ,
              AVERAGE_PERCENTILE( get(paste("ADNI",x,z,sep = "_")) , 
                                  get(paste("bins",x,"gpr",sep = "_")) , 
                                  hv_col_name = "clean_hv_bilateral", hem = "bilateral", xlim = NA ) )



par(mfrow=c(1,1))
boxplot(  c(M_AD , M_MCI , M_CN , F_AD , F_MCI , F_CN ) ~ 
          c( (M_AD * 0 + 1),(M_MCI * 0 + 2), (M_CN * 0 + 3), (F_AD * 0 + 4),(F_MCI * 0 + 5), (F_CN * 0 + 6)  ), 
          col=c( rep(rgb(0,0,1,0.2) , 3), rep(rgb(1,0,0,0.6) , 3) ),
          outline=FALSE , ylab = "Percentile" , xlab = "Strata of Data" , main = "HV Percentile distribution in GPR-Method Nomogram")

abline(h=50 , col='grey' , lty=1)



sep <- SEPERATE( ADNI_male , 'PRS_TH_1' , 0.5)
ADNI_male_low <- sep$Table_low
ADNI_male_high <- sep$Table_high

ADNI_male_high_AD <- filter(ADNI_male_high , VISCODE.y=="scmri" & DX=="Dementia" )
ADNI_male_high_MCI <- filter(ADNI_male_high , VISCODE.y=="scmri" & DX=="MCI" )
ADNI_male_high_CN <- filter(ADNI_male_high , VISCODE.y=="scmri" & DX=="CN" )
ADNI_male_high_NA <- filter(ADNI_male_high , VISCODE.y=="scmri" & is.na(DX) )

ADNI_male_low_AD <- filter(ADNI_male_low , VISCODE.y=="scmri" & DX=="Dementia" )
ADNI_male_low_MCI <- filter(ADNI_male_low , VISCODE.y=="scmri" & DX=="MCI" )
ADNI_male_low_CN <- filter(ADNI_male_low , VISCODE.y=="scmri" & DX=="CN" )
ADNI_male_low_NA <- filter(ADNI_male_low , VISCODE.y=="scmri" & is.na(DX) )


sep <- SEPERATE( ADNI_female , 'PRS_TH_1' , 0.5)
ADNI_female_low <- sep$Table_low
ADNI_female_high <- sep$Table_high

ADNI_female_high_AD <- filter(ADNI_female_high , VISCODE.y=="scmri" & DX=="Dementia" )
ADNI_female_high_MCI <- filter(ADNI_female_high , VISCODE.y=="scmri" & DX=="MCI" )
ADNI_female_high_CN <- filter(ADNI_female_high , VISCODE.y=="scmri" & DX=="CN" )
ADNI_female_high_NA <- filter(ADNI_female_high , VISCODE.y=="scmri" & is.na(DX) )

ADNI_female_low_AD <- filter(ADNI_female_low , VISCODE.y=="scmri" & DX=="Dementia" )
ADNI_female_low_MCI <- filter(ADNI_female_low , VISCODE.y=="scmri" & DX=="MCI" )
ADNI_female_low_CN <- filter(ADNI_female_low , VISCODE.y=="scmri" & DX=="CN" )
ADNI_female_low_NA <- filter(ADNI_female_low , VISCODE.y=="scmri" & is.na(DX) )

# create a bunch of variables with a loop rather than write them all out
# variables are in the form of:
#M_CN_HN <- AVERAGE_PERCENTILE( ADNI_male_high_CN , bins_male , "clean_hv_left" , 1 )
# and we are generating:
#  M_CN_HN, M_CN_H, M_CN_LN, M_CN_L, 
#  F_CN_HN, F_CN_H, F_CN_LN, F_CN_L, 
for( x in c("male","female") )
    for(z in c("AD","MCI","CN","NA") )
      for(j in c("high","low")){
        assign( paste(toupper(substring(x,1,1)),"_",z,"_",toupper(substring(j,1,1)),"N",sep = "") ,
                AVERAGE_PERCENTILE( get(paste("ADNI_",x,"_",j,"_",z,sep = "")) , 
                                    get(paste("bins_",x,"",sep = "")) , 
                                    "clean_hv_bilateral", hem = "bilateral" ) )
        assign( paste(toupper(substring(x,1,1)),"_",z,"_",toupper(substring(j,1,1)),sep = "") ,
                AVERAGE_PERCENTILE( get(paste("ADNI_",x,"_",j,"_",z,sep = "")) , 
                                    get(paste("bins_",x,"_",j,sep = "")) , 
                                    "clean_hv_bilateral", hem = "bilateral" ) )
      }


par(mfrow=c(1,1))

boxplot(  c(M_CN_HN , M_CN_H , M_CN_LN , M_CN_L , F_CN_HN , F_CN_H , F_CN_LN , F_CN_L ) ~ 
            c( (M_CN_HN*0 + 1) , (M_CN_H*0 + 2), (M_CN_LN*0 + 3),(M_CN_L*0 + 4),
               (F_CN_HN*0 + 5) , (F_CN_H*0 + 6),(F_CN_LN*0 + 7),(F_CN_L*0 + 8) ),
          col= c( rep(c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2) ), 2) ),
          outline=FALSE, ylab = "Percentile" , xlab = "Strata of Data" , main = "CN HV Percentile distribution in GPR-Method Nomogram"
)
abline(h=50 , col='grey' , lty=1)

boxplot(  c(M_MCI_HN , M_MCI_H , M_MCI_LN , M_MCI_L ,
            F_MCI_HN , F_MCI_H , F_MCI_LN , F_MCI_L ) ~ 
            c( (M_MCI_HN*0 + 1),(M_MCI_H*0 + 2),(M_MCI_LN*0 + 3),(M_MCI_L*0 + 4),
               (F_MCI_HN*0 + 5),(F_MCI_H*0 + 6),(F_MCI_LN*0 + 7),(F_MCI_L*0 + 8) ) ,
          col= c( rep(c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2) ), 2) ),
          outline=FALSE , ylab = "Percentile" , xlab = "Strata of Data" 
)

boxplot(  c(M_AD_HN , M_AD_H , M_AD_LN , M_AD_L ,
            F_AD_HN , F_AD_H , F_AD_LN , F_AD_L  ) ~ 
            c( (1:length(M_AD_HN) * 0 + 1),(1:length(M_AD_H) * 0 + 2),(1:length(M_AD_LN) * 0 + 3),(1:length(M_AD_L) * 0 + 4),
               (1:length(F_AD_HN) * 0 + 5),(1:length(F_AD_H) * 0 + 6),(1:length(F_AD_LN) * 0 + 7),(1:length(F_AD_L) * 0 + 8) ) ,
          col= c( rep(c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2) ), 2) ),
          outline=FALSE , ylab = "Percentile" , xlab = "Strata of Data" , main = "AD HV Percentile distribution in GPR-Method Nomogram + ADNI CN"
)


par(mfrow=c(1,3))
PLOT_NOMOGRAM_ADNI(bins_male_gpr , ADNI_male , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Male Bilateral" , hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_high , ADNI_male_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs", hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_low , ADNI_male_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "bilateral")

par(mfrow=c(3,6))
PLOT_NOMOGRAM_ADNI(bins_male_gpr , ADNI_male , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Male Left" , hem = "left")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_high , ADNI_male_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs nomogram", hem = "left")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_low , ADNI_male_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "left")
PLOT_NOMOGRAM_ADNI(bins_female_gpr , ADNI_female , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Female Left", hem = "left")
PLOT_NOMOGRAM_ADNI(bins_female_gpr_high , ADNI_female_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs", hem = "left")
PLOT_NOMOGRAM_ADNI(bins_female_gpr_low , ADNI_female_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "left")

PLOT_NOMOGRAM_ADNI(bins_male_gpr , ADNI_male , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Male Right" , hem = "right")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_high , ADNI_male_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs", hem = "right")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_low , ADNI_male_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "right")
PLOT_NOMOGRAM_ADNI(bins_female_gpr , ADNI_female , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Female Right", hem = "right")
PLOT_NOMOGRAM_ADNI(bins_female_gpr_high , ADNI_female_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs", hem = "right")
PLOT_NOMOGRAM_ADNI(bins_female_gpr_low , ADNI_female_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "right")

PLOT_NOMOGRAM_ADNI(bins_male_gpr , ADNI_male , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Male Bilateral" , hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_high , ADNI_male_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs", hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_male_gpr_low , ADNI_male_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_female_gpr , ADNI_female , ylim = c(1600,5200) , xlim=c(45,90) , title = "Full Nomogram Female Bilateral", hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_female_gpr_high , ADNI_female_high , ylim = c(1600,5200) , xlim=c(45,90), title = "high prs", hem = "bilateral")
PLOT_NOMOGRAM_ADNI(bins_female_gpr_low , ADNI_female_low , ylim = c(1600,5200) , xlim=c(45,90) , title = "low prs", hem = "bilateral")


##### LONGITUDINAL ANALYSIS #########

par(mfrow=c(1,3))
IDS <- c(4432 , 4380 , 4960 , 4042 , 4631 , 4513 , 2274 , 2195 , 4406 , 5135 , 4706 , 2216 , 4188 , 2378 , 2195 , 4885 , 4218)
PLOT_NOMOGRAM_ADNI_LONGITUDINAL( bins_male_gpr , ADNI_male , hem="bilateral" , title=paste("Male HV"), xlim = c(45,100), ylim=c(1500,5200) , IDS = NA)
PLOT_NOMOGRAM_ADNI_LONGITUDINAL( bins_male_gpr_high , ADNI_male_high , hem="bilateral" , title=paste("High PRS"), xlim = c(45,100), ylim=c(1500,5200) ,IDS = NA)
PLOT_NOMOGRAM_ADNI_LONGITUDINAL( bins_male_gpr_low , ADNI_male_low , hem="bilateral" , title=paste("Low PRS"), xlim = c(45,100), ylim=c(1500,5200) , IDS = NA)

(mean(F_CN_HN, na.rm = TRUE) - mean(F_CN_LN, na.rm = TRUE) +
  mean(M_CN_HN, na.rm = TRUE) - mean(M_CN_LN, na.rm = TRUE)  ) / 2


(mean(F_CN_H, na.rm = TRUE) - mean(F_CN_L, na.rm = TRUE) +
  mean(M_CN_H, na.rm = TRUE) - mean(M_CN_L, na.rm = TRUE) ) / 2

# code used to generate figure 1 from AAIC abstract
par(mfrow=c(1,4))
PLOT_NOMOGRAM( bins_male , hem="bilateral" , paste("SWM Nomogram") , xlim=c(52,78) , ylim=c(3100 , 5000) )
PLOT_NOMOGRAM( bins_male_gpr_s_smth , hem="bilateral" , title = "GPR Nomogram" , xlim=c(52,78) , ylim=c(3100 , 5000) )
PLOT_NOMOGRAM_COMPARE( bins_male_gpr_s_smth , bins_male_gpr_high, hem="bilateral" , title = "GPR Nomogram - High PRS" , xlim=c(52,78) , ylim=c(3100 , 5000))
PLOT_NOMOGRAM_COMPARE( bins_male_gpr_s_smth , bins_male_gpr_low, hem="bilateral" , title = "GPR Nomogram - Low PRS" , xlim=c(52,78) , ylim=c(3100 , 5000))

# code used to generate figure 2 from AAIC abstract
boxplot(  c(M_CN_HN , M_CN_LN , M_CN_H , M_CN_L ) ~ 
            c( (M_CN_HN*0 + 1) , (M_CN_LN*0 + 2), (M_CN_H*0 + 3),(M_CN_L*0 + 4)),
          col= c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2)),
          ylab = "Percentile", xlab="",main = "MALE",
          names=c("High Samples\nRegular Nomogram" , "Low Samples\nRegular Nomogram" , "High Samples\nHigh Nomogram" , "Low Samples\nLow Nomogram"),
          outline = FALSE
)

abline(h=50,col='grey',lty=2)

boxplot(  c(M_AD_HN , M_AD_LN , M_AD_H , M_AD_L ) ~ 
            c( (M_AD_HN*0 + 1) , (M_AD_LN*0 + 2), (M_AD_H*0 + 3),(M_AD_L*0 + 4)),
          col= c( rep(rgb(0,0,1,0.6) , 2), rep(rgb(0,0,1,0.2) , 2)),
          ylab = "Percentile", xlab="",main = "MALE",
          names=c("High Samples\nRegular Nomogram" , "Low Samples\nRegular Nomogram" , "High Samples\nHigh Nomogram" , "Low Samples\nLow Nomogram"),
          outline = FALSE
)


# Exploring impact on MCI

par(mfrow=c(4,4))

MCI_Analysis(ADNI_male_high , "bilateral" , bins_male_gpr )
MCI_Analysis(ADNI_male_low , "bilateral" , bins_male_gpr )
MCI_Analysis(ADNI_male_high , "bilateral" , bins_male_gpr_high )
MCI_Analysis(ADNI_male_low , "bilateral" , bins_male_gpr_low)

MCI_Analysis(ADNI_female_high , "bilateral" , bins_female_gpr )
MCI_Analysis(ADNI_female_low , "bilateral" , bins_female_gpr )
MCI_Analysis(ADNI_female_high , "bilateral" , bins_female_gpr_high )
MCI_Analysis(ADNI_female_low , "bilateral" , bins_female_gpr_low)

##

MCI_Full_1 <- MCI_Analysis_AM(ADNI_female_high , "bilateral" , bins_female_gpr)
MCI_Full_2 <- MCI_Analysis_AM(ADNI_female_low , "bilateral" , bins_female_gpr)
MCI_Full_3 <- MCI_Analysis_AM(ADNI_male_high , "bilateral" , bins_male_gpr )
MCI_Full_4 <- MCI_Analysis_AM(ADNI_male_low , "bilateral" , bins_male_gpr )
MCI_Full_1$sex <- 0
MCI_Full_2$sex <- 0
MCI_Full_3$sex <- 1
MCI_Full_4$sex <- 1

MCI_Full <- rbind( MCI_Full_1 , MCI_Full_2 , MCI_Full_3 , MCI_Full_4)
#summary(glm( convert ~ percent + age + sex , data = MCI_Full , family=binomial) )
summary( coxph( Surv(covert_time , convert) ~ percent + sex + age , data = MCI_Full) )

MCI_PRS_1 <- MCI_Analysis_AM(ADNI_female_high , "bilateral" , bins_female_gpr_high)
MCI_PRS_2 <- MCI_Analysis_AM(ADNI_female_low , "bilateral" , bins_female_gpr_low)
MCI_PRS_3 <- MCI_Analysis_AM(ADNI_male_high , "bilateral" , bins_male_gpr_high)
MCI_PRS_4 <- MCI_Analysis_AM(ADNI_male_low , "bilateral" , bins_male_gpr_low) 
MCI_PRS_1$sex <- 0
MCI_PRS_2$sex <- 0
MCI_PRS_3$sex <- 1
MCI_PRS_4$sex <- 1

MCI_PRS <- rbind( MCI_PRS_1 , MCI_PRS_2 , MCI_PRS_3 , MCI_PRS_4)
#summary(glm( convert ~ percent + age + sex , data = MCI_PRS , family=binomial) )
summary( coxph( Surv(covert_time , convert) ~ percent + sex + age , data = MCI_PRS) )


IDS <- unique(ADNI_table$RID.x)
ADNI_table[ADNI_table$RID.x == IDS , "RID.x"]


plot(ukb_male$AGE_Latest , ukb_male$clean_hv_left , xlim = c(40,90) , ylim=c(1000,7000))
points(ukb_female$AGE_Latest , ukb_female$clean_hv_left )
points(ADNI_male$TRUE.AGE , ADNI_male$clean_hv_left , col = "red" )
points(ADNI_female$TRUE.AGE , ADNI_female$clean_hv_left , col = "red" )


par(mfrow=c(1,1))
ukb_male <- ukb_male[ order(ukb_male[,"AGE_Latest"]), ]

length <- 5
bins_1_gpr <- bins_male_gpr
bins_2_gpr <- GPR_ANALYSIS(ukb_male[floor(seq(1 , nrow(ukb_male), (nrow(ukb_male)/(length-1))-1 )),])

bins_1_swm <- bins_male
bins_2_swm <- WINDOW_ANALYSIS(ukb_male[floor(seq(1 , nrow(ukb_male), (nrow(ukb_male)/(length-1))-1 )),])

# Checking the data effeciency of building nomograms with both methods.
# we want to make a statement like: genrating nomograms with only ()% of the samples produces nomograms that are ()% similar.
res <- c()
res2 <- c()
num_iterations <- 10 

x_vals <- c(seq(2000,5000,1000))
x_vals <- c(seq(10,90,10) , seq(100,1000,100) ) # , seq(1000,2000,1000) )

for( length in x_vals){
  for( i in 1:num_iterations ){
  indices <- floor(seq(1 , nrow(ukb_male), (nrow(ukb_male)/(length-1))-1 ))
  indices[2:(length(indices)-1)] <- floor(jitter(indices[2:(length(indices)-1)], amount = 3))
  samples <- ukb_male[indices,]
  
  bins_2_gpr <- GPR_ANALYSIS( samples , XX = ages_male)
  bins_2_swm <- WINDOW_ANALYSIS(samples)
  
  res2 <- c(res2 , ( (NOMOGRAM_DIFF_INTERPOLATE(bins_1_gpr , bins_2_gpr , age_range = range(samples$AGE_Latest) ) ) ) )
  res <- c(res , ( (NOMOGRAM_DIFF_INTERPOLATE(bins_1_swm , bins_2_swm , age_range = range(samples$AGE_Latest) ) ) ) ) 
  }
}

means <- apply(matrix(res, nrow=num_iterations,ncol=length(x_vals)) , 2 , mean)
sds <- apply(matrix(res, nrow=num_iterations,ncol=length(x_vals)) , 2 , sd)

means2 <- apply(matrix(res2, nrow=num_iterations,ncol=length(x_vals)) , 2 , mean)
sds2 <- apply(matrix(res2, nrow=num_iterations,ncol=length(x_vals)) , 2 , sd)

plot(x_vals , means , ylim = c(0,300) , col="blue" , type = "l")
points(x_vals , means2 , ylim = c(0,300) , col="red" , type = "l")

lines( x_vals , ( means +  qnorm(0.975) * sds / 10 ) )
lines( x_vals , ( means -  qnorm(0.975) * sds / 10 ) )

lines( x_vals , ( means2 +  qnorm(0.975) * sds2 / 10 ) )
lines( x_vals , ( means2 -  qnorm(0.975) * sds2 / 10 ) )



times <- c()
x_vals <- rep(900,100)
for( length in x_vals){
    indices <- floor(seq(1 , nrow(ukb_male), (nrow(ukb_male)/(length-1))-1 ))
    samples <- ukb_male[indices,]
    start_time <- Sys.time()
    bins_2_gpr <- GPR_ANALYSIS( samples , XX = ages_male)
    end_time <- Sys.time()
    duration <- as.numeric( end_time - start_time , units="secs")
    times <- c(times , duration)
}



# plotting the PRS distribution in both datasets to check if harmonization between them is needed
plot( -60:60 , dnorm(-60:60 , mean = mean(ukb_male$PRS_TH_1,na.rm = TRUE) , sd = sd(ukb_male$PRS_TH_1,na.rm = TRUE) ) , type = "l" , col="red" , ylab = "" , xlab = "PRS value" , main = "Gaussian Fitting of PRS scores")
lines( -60:60 , dnorm(-60:60 , mean = mean(ADNI_male$PRS_TH_1,na.rm = TRUE) , sd = sd(ADNI_male$PRS_TH_1,na.rm = TRUE) ) , type = "l" , col="blue")
legend(1,95,legend = c("UKBB","ADNI") , col = c("red","blue"), lty=c(1,1))
legend(1,95,legend = c("UKBB","ADNI") , col = c("red","blue"), lty=c(1,1) , cex=0.8)
legend(-30,0.15,legend = c("UKBB","ADNI") , col = c("red","blue"), lty=c(1,1) , cex=0.8)




par(mfrow=c(2,3))

indices <- floor(seq(1 , nrow(ukb_male), nrow(ukb_male)/20000 ))
ukb_male <- ukb_male[ order(ukb_male[,"AGE_Latest"]), ]
full_table <- ukb_male[indices,]

#indices <- floor(seq(1 , nrow(ukb_female), nrow(ukb_female)/20000 ))
#ukb_female <- ukb_female[ order(ukb_female[,"AGE_Latest"]), ]
#full_table <- ukb_female[indices,]


aged_table <- full_table[ (full_table$AGE_Latest >= 46 & full_table$AGE_Latest <= 66), ]
table_1 <- aged_table[ sample( 1:nrow(aged_table), 1000 , replace = FALSE), ]
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

diffs_1 <- c()
diffs_2 <- c()
diffs_3 <- c()
for( i in 30:90 ){ 
  diffs_1 <- c( diffs_1 , NOMOGRAM_DIFF_INTERPOLATE(bins_1,bins_female_expanded_age, age_range = c(i,i+1)) )
  diffs_2 <- c( diffs_2 , NOMOGRAM_DIFF_INTERPOLATE(bins_2,bins_female_expanded_age, age_range = c(i,i+1)) )
  diffs_3 <- c( diffs_3 , NOMOGRAM_DIFF_INTERPOLATE(bins_3,bins_female_expanded_age, age_range = c(i,i+1)) )
}

plot( 30:90 , diffs_1 , xlab = "Age" , ylab = "Volume Difference" , main = "46-66 range")
abline(v=46)
abline(v=66)

plot( 30:90 , diffs_2 , xlab = "Age" , ylab = "Volume Difference" , main = "50-70 range")
abline(v=50)
abline(v=70)

plot( 30:90 , diffs_3 , xlab = "Age" , ylab = "Volume Difference" , main = "60-80 range")
abline(v=60)
abline(v=80)

mean( mean(diffs_1[16:36]) , mean(diffs_2[20:40]) , mean(diffs_3[30:50]) )

S1 <- rev((diffs_1[1:16] + diffs_2[5:20] + diffs_3[15:30]) / 3 )
S2 <- (diffs_3[50:61] + diffs_2[40:51] + diffs_1[36:47] ) / 3

SS <- (head(S1, length(S2)) + S2)/2

mean(SS)





aged_table <- ukb_male[ (ukb_male$AGE_Latest >= 50 & ukb_male$AGE_Latest <= 80), ]

table_1 <- aged_table[seq(1,nrow(aged_table),nrow(aged_table)/600),]
PLOT_NOMOGRAM( bins_1 <- GPR_ANALYSIS(table_1) , title = "600 samples")
abline(v=max(table_1$AGE_Latest) , col="black")
abline(v=min(table_1$AGE_Latest) , col="black")

table_2 <- aged_table[seq(1,nrow(aged_table),nrow(aged_table)/1000),]
PLOT_NOMOGRAM( bins_2 <- GPR_ANALYSIS(table_2) , title = "1000 samples")
abline(v=max(table_2$AGE_Latest) , col="black")
abline(v=min(table_2$AGE_Latest) , col="black")

table_3 <- aged_table[seq(1,nrow(aged_table),nrow(aged_table)/2000),]
PLOT_NOMOGRAM( bins_3 <- GPR_ANALYSIS(table_3) , title= "2000 samples")
abline(v=max(table_3$AGE_Latest) , col="black")
abline(v=min(table_3$AGE_Latest) , col="black")

table_4 <- aged_table[seq(1,nrow(aged_table),nrow(aged_table)/5000),]
PLOT_NOMOGRAM( bins_4 <- GPR_ANALYSIS(table_4) , title = "5000 samples")
abline(v=max(table_4$AGE_Latest) , col="black")
abline(v=min(table_4$AGE_Latest) , col="black")

diffs_1 <- c()
diffs_2 <- c()
diffs_3 <- c()
diffs_4 <- c()
for( i in 30:90 ){ 
  diffs_1 <- c( diffs_1 , NOMOGRAM_DIFF_INTERPOLATE(bins_1,bins_male_expanded_age, age_range = c(i,i+1)) )
  diffs_2 <- c( diffs_2 , NOMOGRAM_DIFF_INTERPOLATE(bins_2,bins_male_expanded_age, age_range = c(i,i+1)) )
  diffs_3 <- c( diffs_3 , NOMOGRAM_DIFF_INTERPOLATE(bins_3,bins_male_expanded_age, age_range = c(i,i+1)) )
  diffs_4 <- c( diffs_4 , NOMOGRAM_DIFF_INTERPOLATE(bins_4,bins_male_expanded_age, age_range = c(i,i+1)) )
}

plot( 30:90 , diffs_1 , xlab = "Age" , ylab = "Volume Difference" , ylim = c(0,16))
abline(v=max(table_1$AGE_Latest) , col="black")
abline(v=min(table_1$AGE_Latest) , col="black")

plot( 30:90 , diffs_2 , xlab = "Age" , ylab = "Volume Difference" , ylim = c(0,16))
abline(v=max(table_2$AGE_Latest) , col="black")
abline(v=min(table_2$AGE_Latest) , col="black")

plot( 30:90 , diffs_3 , xlab = "Age" , ylab = "Volume Difference", ylim = c(0,16))
abline(v=max(table_3$AGE_Latest) , col="black")
abline(v=min(table_3$AGE_Latest) , col="black")

plot( 30:90 , diffs_4 , xlab = "Age" , ylab = "Volume Difference", ylim = c(0,16))
abline(v=max(table_4$AGE_Latest) , col="black")
abline(v=min(table_4$AGE_Latest) , col="black")

mean(diffs_1[20:50])
mean(diffs_2[20:50])
mean(diffs_3[20:50])
mean(diffs_4[20:50])

diffs_1[50:60]




thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.06','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')

par(mfrow=c(3,5))
for(TH in thresholds){
  hist(ukb_male[,TH] , col = alpha("blue",0.3) , xlim = range(ADNI_male[,TH]) , n=50 , main = TH )
  hist(ADNI_male[,TH] , col = alpha("red",0.3) , xlim = range(ADNI_male[,TH]) , n=100 , add=TRUE)
}


par(mfrow=c(2,3))
PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_high , xlim=c(55,75) , ylim=c(3200,5100) , title = "HV HIGH")
PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_high_ICV , xlim=c(55,75), ylim=c(3200,5100) , title = "ICV HIGH")
PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_high_AD , xlim=c(55,75), ylim=c(3200,5100) , title = "AD HIGH")

PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_low , xlim=c(55,75), ylim=c(3200,5100) , title = "HV LOW")
PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_low_ICV , xlim=c(55,75), ylim=c(3200,5100) , title = "ICV LOW")
PLOT_NOMOGRAM_COMPARE(bins_male , bins_male_low_AD , xlim=c(55,75), ylim=c(3200,5100) , title = "AD LOW")
