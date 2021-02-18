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

bins_male_gpr <- GPR_ANALYSIS( ukb_male[1:1000 , ] )
bins_female_gpr <- GPR_ANALYSIS( ukb_female[1:1000 , ] )

######### PLOTTING  NOMOGRAMS ######### 

par(mfrow=c(2,2))

PLOT_NOMOGRAM( bins_male , "bilateral" , "SW male HV")
PLOT_NOMOGRAM( bins_female , "bilateral" , "SW Female HV")

PLOT_NOMOGRAM( bins_male_gpr , "bilateral" , "GPR male Left HV")
PLOT_NOMOGRAM( bins_female_gpr , "bilateral" , "GPR Female Left HV")


###### PRS ANALYSIS #######

MAKE_PRS_PLOTS( ukb_male , ukb_female , 
                confounder_columns = c("AGE_Latest","Sex","Genetic.PC.1","Genetic.PC.2","Genetic.PC.3") )


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

bins_male_gpr_low  <- GPR_ANALYSIS( ukb_male_low[1:1000,] )
bins_male_gpr_high <- GPR_ANALYSIS( ukb_male_high[1:1000,] ) 

bins_female_gpr_low  <- GPR_ANALYSIS( ukb_female_low[1:1000,] )
bins_female_gpr_high <- GPR_ANALYSIS( ukb_female_high[1:1000,] ) 

# plotting is based on the variables above. to avoid writing many lines, 
# I loop through variable names, and get the value of the variable in each loop

# one plot to show all separate nomograms 
par(mfrow=c(2,3))
for( sex in c("male" , "female"))
  for( strat in c("" , "_low" , "_high") ){
      bins <- get(paste("bins_",sex,strat,sep=""))
      title <- paste("HV",toupper(sex),strat,"prs")
      PLOT_NOMOGRAM( bins , hem="bilateral" , title=title )
  }

# one more plot with nomograms overlaid on each other 
par(mfrow=c(1,3))
for( sex in c("male" , "female")){
  bins_1 <- get(paste("bins",sex,sep="_"))
  PLOT_NOMOGRAM( bins_1 , xlim=c(51,80) , title=paste(toupper(sex),"NOMOGRAM"))
  for( strat in c("high" , "low") ){
      bins_2 <- get(paste("bins",sex,strat,sep="_"))
      title <- paste(toupper(sex),"Normal vs",strat,"prs")
      PLOT_NOMOGRAM_COMPARE( bins_1 , bins_2 , hem = "bilateral" , "black" , "red" , title, xlim=c(51,80))
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
PLOT_NOMOGRAM_ADNI(bins_male_gpr , ADNI_male , hem="bilateral" , title=paste("Nomogram Male HV") , xlim = c(52,80) , ylim=c(2000,5200) )
PLOT_NOMOGRAM_ADNI(bins_female_gpr , ADNI_female , hem="bilateral" , title=paste("Nomogram Female HV") , xlim = c(52,80) , ylim=c(2000,5200) )



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
                                    get(paste("bins_",x,"_gpr",sep = "")) , 
                                    "clean_hv_bilateral", hem = "bilateral" ) )
        assign( paste(toupper(substring(x,1,1)),"_",z,"_",toupper(substring(j,1,1)),sep = "") ,
                AVERAGE_PERCENTILE( get(paste("ADNI_",x,"_",j,"_",z,sep = "")) , 
                                    get(paste("bins_",x,"_gpr_",j,sep = "")) , 
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

##### LONGITUDINAL ANALYSIS #########

par(mfrow=c(1,3))

PLOT_NOMOGRAM_ADNI_LONGITUDINAL( bins_female_gpr , ADNI_female , hem="bilateral" , title=paste("Female HV"), xlim = c(52,100), ylim=c(2000,5200) )
PLOT_NOMOGRAM_ADNI_LONGITUDINAL( bins_female_gpr_high , ADNI_female_high , hem="bilateral" , title=paste("Female HV - High PRS"), xlim = c(52,80), ylim=c(2000,5200) )
PLOT_NOMOGRAM_ADNI_LONGITUDINAL( bins_female_gpr_low , ADNI_female_low , hem="bilateral" , title=paste("Female HV - Low PRS"), xlim = c(52,80), ylim=c(2000,5200) )

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

MCI_Analysis(ADNI_male_high , "bilateral" , bins_male_gpr , "Male High-Full")
MCI_Analysis(ADNI_male_low , "bilateral" , bins_male_gpr , "Male Low-Full")
MCI_Analysis(ADNI_male_high , "bilateral" , bins_male_gpr_high , "Male High-High")
MCI_Analysis(ADNI_male_low , "bilateral" , bins_male_gpr_low , "Male Low-Low")

MCI_Analysis(ADNI_female_high , "bilateral" , bins_female_gpr , "Female High-Full")
MCI_Analysis(ADNI_female_low , "bilateral" , bins_female_gpr , "Female Low-Full")
MCI_Analysis(ADNI_female_high , "bilateral" , bins_female_gpr_high , "Female High-High")
MCI_Analysis(ADNI_female_low , "bilateral" , bins_female_gpr_low , "Female Low-Low")

##

MCI_Full_1 <- MCI_Analysis(ADNI_female_high , "bilateral" , bins_female_gpr, "Female High-Full")
MCI_Full_2 <- MCI_Analysis(ADNI_female_low , "bilateral" , bins_female_gpr , "Female Low-Full")
MCI_Full_3 <- MCI_Analysis(ADNI_male_high , "bilateral" , bins_male_gpr , "Male High-Full")
MCI_Full_4 <- MCI_Analysis(ADNI_male_low , "bilateral" , bins_male_gpr , "Male Low-Full")
MCI_Full_1$sex <- 0
MCI_Full_2$sex <- 0
MCI_Full_3$sex <- 1
MCI_Full_4$sex <- 1

MCI_Full <- rbind( MCI_Full_1 , MCI_Full_2 , MCI_Full_3 , MCI_Full_4)
summary(glm( convert ~ percent + age + sex , data = MCI_Full , family=binomial) )

MCI_PRS_1 <- MCI_Analysis(ADNI_female_high , "bilateral" , bins_female_gpr_high , "Female High-High")
MCI_PRS_2 <- MCI_Analysis(ADNI_female_low , "bilateral" , bins_female_gpr_low , "Female Low-Low")
MCI_PRS_3 <- MCI_Analysis(ADNI_male_high , "bilateral" , bins_male_gpr_high , "Male High-High")
MCI_PRS_4 <- MCI_Analysis(ADNI_male_low , "bilateral" , bins_male_gpr_low , "Male Low-Low") 
MCI_PRS_1$sex <- 0
MCI_PRS_2$sex <- 0
MCI_PRS_3$sex <- 1
MCI_PRS_4$sex <- 1

MCI_PRS <- rbind( MCI_PRS_1 , MCI_PRS_2 , MCI_PRS_3 , MCI_PRS_4)
summary(glm( convert ~ percent + age + sex , data = MCI_PRS , family=binomial) )


IDS <- unique(ADNI_table$RID.x)
ADNI_table[ADNI_table$RID.x == IDS , "RID.x"]

