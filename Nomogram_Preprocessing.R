# Script to prep UK Biobank Data for Nomogram generation
# By: Mohammed Janahi


##### SUPPORT FUNCTIONS  ########### 

# Function to perform the sliding window analysis
# 
# INPUTS:
# 
#
# OUTPUT:
#
#
# ASSUMPTIONS:
#
#
WINDOW_ANALYSIS <- function( ukb , percent_samples_per_bin=10 ){
  
  NUM_SAMPLES <- nrow(ukb)
  samples_per_bin <- NUM_SAMPLES / percent_samples_per_bin
  
  one_percent_of_samples <- NUM_SAMPLES/100
  
  bins <- data.frame(min_age=double(), max_age=double(), mid_age=double(), mean=double(), std=double(),
                     q2.5=double(), q5=double(), q10=double(), q25=double(), q50=double(), q75=double(), q90=double(), q95=double(), q97.5=double(), 
                     mean_left=double(), std_left=double(), 
                     l.q2.5=double(), l.q5=double(), l.q10=double(), l.q25=double(), l.q50=double(), l.q75=double(), l.q90=double(), l.q95=double(), l.q97.5=double(), 
                     mean_right=double(), std_right=double(),
                     r.q2.5=double(), r.q5=double(), r.q10=double(), r.q25=double(), r.q50=double(), r.q75=double(), r.q90=double(), r.q95=double(), r.q97.5=double()
                     )
  
  for (i in seq( 1 , (NUM_SAMPLES - samples_per_bin) , one_percent_of_samples) ){
    
    #window_of_samples_left <- ukb[i:(i+samples_per_bin) , c("X21003.2.0","resid_hv_volumes_left")]
    #window_of_samples_right <- ukb[i:(i+samples_per_bin) , c("X21003.2.0","resid_hv_volumes_right")]
    hv_left_col_name  <- "clean_hv_volumes_left"  #"X25019.2.0"
    hv_right_col_name <- "clean_hv_volumes_right" #"X25020.2.0"
    
    window_of_samples_left <- ukb[i:(i+samples_per_bin) , c("Age.Imaging.Visit" , "HV.Left" , hv_left_col_name)]
    window_of_samples_right <- ukb[i:(i+samples_per_bin) , c("Age.Imaging.Visit" , "HV.Right" , hv_right_col_name)]
    
    win_ages      <- window_of_samples_left$Age.Imaging.Visit
    win_min_age   <- min(win_ages  , na.rm = TRUE)
    win_max_age   <- max(win_ages  , na.rm = TRUE)
    win_mean_age  <- mean(win_ages , na.rm = TRUE)
    
    win_hv_left   <- window_of_samples_left$HV.Left
    win_hv_right  <- window_of_samples_right$HV.Right
    win_hv_both   <- unlist(c( win_hv_left , win_hv_right ))
    
    win_hv_left   <- window_of_samples_left[,hv_left_col_name]
    win_hv_right  <- window_of_samples_right[,hv_right_col_name]
    win_hv_both   <- unlist(c( win_hv_left , win_hv_right ))
    
    win_mean_hv   <- mean(win_hv_both, na.rm = TRUE) 
    win_std_hv    <- sd(win_hv_both , na.rm = TRUE)
    win_quantiles <- quantile( (win_hv_both) , probs = c(0.025, 0.05 , 0.1 , 0.25 , 0.5 , 0.75 , 0.9 , 0.95, 0.975) , na.rm = TRUE) 
                     # + win_mean_hv
    
    win_mean_hv_left   <- mean(win_hv_left, na.rm = TRUE) 
    win_std_hv_left    <- sd(win_hv_left , na.rm = TRUE)
    win_quantiles_left <- quantile( (win_hv_left) , probs = c(0.025, 0.05 , 0.1 , 0.25 , 0.5 , 0.75 , 0.9 , 0.95, 0.975) , na.rm = TRUE) 
                          # + win_mean_hv_left
    
    win_mean_hv_right  <- mean(win_hv_right, na.rm = TRUE) 
    win_std_hv_right   <- sd(win_hv_right , na.rm = TRUE)
    win_quantiles_right <- quantile( (win_hv_right) , probs = c(0.025, 0.05 , 0.1 , 0.25 , 0.5 , 0.75 , 0.9 , 0.95, 0.975) , na.rm = TRUE) 
                           # + win_mean_hv_right
    

    row <- unlist(list( win_min_age , win_max_age , win_mean_age , win_mean_hv , win_std_hv , win_quantiles ,
                        win_mean_hv_left , win_std_hv_left , win_quantiles_left , 
                        win_mean_hv_right , win_std_hv_right, win_quantiles_right ))
    
    bins[nrow(bins) + 1,]  = row
  }
  return(bins)
}

# Function to plot nomograms given sliding window bins object
#
PLOT_NOMOGRAM <- function( bins , win , left , title ){
  
  centers <- bins$mid_age
  mn <- (bins$mean_left * left )  + (bins$mean_right * abs(left-1) ) # Instead of an if statement, weight by 0/1 for correct side
  q2.5 <- (bins$l.q2.5 * left )  + (bins$r.q2.5 * abs(left-1) ) 
  q5   <- (bins$l.q5 * left )    + (bins$r.q5 * abs(left-1) ) 
  q10  <- (bins$l.q10 * left )   + (bins$r.q10 * abs(left-1) )
  q25  <- (bins$l.q25 * left )   + (bins$r.q25 * abs(left-1) )
  q50  <- (bins$l.q50 * left )   + (bins$r.q50 * abs(left-1) )
  q75  <- (bins$l.q75 * left )   + (bins$r.q75 * abs(left-1) )
  q90  <- (bins$l.q90 * left )   + (bins$r.q90 * abs(left-1) )
  q95  <- (bins$l.q95 * left )   + (bins$r.q95 * abs(left-1) )
  q97.5<- (bins$l.q97.5 * left ) + (bins$r.q97.5 * abs(left-1) )
  
  smoothed_mn <- smth( mn, method = 'gaussian', window=win)
  
  min_x <- floor(centers[ 1 + floor(win/2) ])
  max_x <- floor(centers[ length(centers) - floor(win/2)])
  
  plot( main = paste(title," smoothing kernel width = ",win) , xlab = "Age" , ylab = "HV" ,
        centers , smoothed_mn,
        ylim=c(2400,5000), xlim=c(  min_x , max_x ) , 
        xaxt="n" , yaxt="none", type='l')
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(2400,5000,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(2400, 5000, 100) , col='grey' , lty=2)
  
  lines( centers , smth( q2.5,  method = 'gaussian', window=win) )
  lines( centers , smth( q5,    method = 'gaussian', window=win) )
  lines( centers , smth( q10,   method = 'gaussian', window=win) )
  lines( centers , smth( q25,   method = 'gaussian', window=win) )
  lines( centers , smth( q75,   method = 'gaussian', window=win) )
  lines( centers , smth( q90,   method = 'gaussian', window=win) )
  lines( centers , smth( q95,   method = 'gaussian', window=win) )
  lines( centers , smth( q97.5, method = 'gaussian', window=win) )
  
  str_ind <- (1 + floor(win/2))
  end_ind <- (length(centers) - floor(win/2))
  centers_sml = centers[str_ind:end_ind]
  q2.5_sml <- smth( q2.5,   method = 'gaussian', window=win)[str_ind:end_ind]
  q97.5_sml <- smth( q97.5,   method = 'gaussian', window=win)[str_ind:end_ind]
  
  q5_sml <- smth( q5,   method = 'gaussian', window=win)[str_ind:end_ind]
  q95_sml <- smth( q95,   method = 'gaussian', window=win)[str_ind:end_ind]
  
  q10_sml <- smth( q10,   method = 'gaussian', window=win)[str_ind:end_ind]
  q90_sml <- smth( q90,   method = 'gaussian', window=win)[str_ind:end_ind]
  
  q25_sml <- smth( q25,   method = 'gaussian', window=win)[str_ind:end_ind]
  q75_sml <- smth( q75,   method = 'gaussian', window=win)[str_ind:end_ind]
  
  polygon( c( centers_sml , rev(centers_sml) ) , c(q2.5_sml , rev(q97.5_sml)) , col=rgb(0,1,0,0.1) , border = NA)
  polygon( c( centers_sml , rev(centers_sml) ) , c(q5_sml , rev(q95_sml)) , col=rgb(0,1,0,0.1) , border = NA)
  polygon( c( centers_sml , rev(centers_sml) ) , c(q10_sml , rev(q90_sml)) , col=rgb(0,1,0,0.1) , border = NA)
  polygon( c( centers_sml , rev(centers_sml) ) , c(q25_sml , rev(q75_sml)) , col=rgb(0,1,0,0.1) , border = NA)
  
}

###### MAIN CODE STARTS HERE #######

library(smoother)
library(TTR)

# UKB Table is saved as CSV file on lacal machine
# Two files to load because we selected the columns in two batches
ukb <- read.csv("~/Desktop/ukb_data/ukb_table.csv")
ukb_ext <- read.csv("~/Desktop/ukb_data/ukb_features_v2.csv")
# Merge the two tables
ukb <- merge(ukb , ukb_ext , by.x="eid" , by.y="eid")

#ukb <- read.csv("~/Desktop/ukb_data/with_imaging_40k/ukb_features_042020.csv")

# For readability, rename some columns
names(ukb)[names(ukb) == 'X31.0.0']    <- 'Sex'
names(ukb)[names(ukb) == 'X21003.0.0'] <- 'Age.When.Attended.Assesment.Center'
names(ukb)[names(ukb) == 'X21003.2.0'] <- 'Age.Imaging.Visit'
names(ukb)[names(ukb) == 'X25000.2.0'] <- 'Volumetric.scaling.from.T1.head.image.to.standard.space'
names(ukb)[names(ukb) == 'X53.2.0']    <- 'Scan.Date'
names(ukb)[names(ukb) == 'X25019.2.0'] <- 'HV.Left'
names(ukb)[names(ukb) == 'X25020.2.0'] <- 'HV.Right'
names(ukb)[names(ukb) == 'X21000.0.0'] <- 'Ethnic.Background'
names(ukb)[names(ukb) == 'X22009.0.1'] <- 'Genetic.PC.1'
names(ukb)[names(ukb) == 'X22009.0.2'] <- 'Genetic.PC.2'
names(ukb)[names(ukb) == 'X22009.0.3'] <- 'Genetic.PC.3'
names(ukb)[names(ukb) == 'X22009.0.4'] <- 'Genetic.PC.4'
names(ukb)[names(ukb) == 'X22009.0.5'] <- 'Genetic.PC.5'

names(ukb)[names(ukb) == 'X20002.2.0']  <- 'Self.Reported.Conditions.1'
names(ukb)[names(ukb) == 'X20002.2.1']  <- 'Self.Reported.Conditions.2'
names(ukb)[names(ukb) == 'X20002.2.2']  <- 'Self.Reported.Conditions.3'
names(ukb)[names(ukb) == 'X20002.2.3']  <- 'Self.Reported.Conditions.4'
names(ukb)[names(ukb) == 'X20002.2.4']  <- 'Self.Reported.Conditions.5'
names(ukb)[names(ukb) == 'X20002.2.5']  <- 'Self.Reported.Conditions.6'
names(ukb)[names(ukb) == 'X20002.2.6']  <- 'Self.Reported.Conditions.7'
names(ukb)[names(ukb) == 'X20002.2.7']  <- 'Self.Reported.Conditions.8'
names(ukb)[names(ukb) == 'X20002.2.8']  <- 'Self.Reported.Conditions.9'
names(ukb)[names(ukb) == 'X20002.2.9']  <- 'Self.Reported.Conditions.10'
names(ukb)[names(ukb) == 'X20002.2.10'] <- 'Self.Reported.Conditions.11'
names(ukb)[names(ukb) == 'X20002.2.11'] <- 'Self.Reported.Conditions.12'
names(ukb)[names(ukb) == 'X20002.2.12'] <- 'Self.Reported.Conditions.13'
names(ukb)[names(ukb) == 'X20002.2.13'] <- 'Self.Reported.Conditions.14'
names(ukb)[names(ukb) == 'X20002.2.14'] <- 'Self.Reported.Conditions.15'
names(ukb)[names(ukb) == 'X20002.2.15'] <- 'Self.Reported.Conditions.16'
names(ukb)[names(ukb) == 'X20002.2.16'] <- 'Self.Reported.Conditions.17'
names(ukb)[names(ukb) == 'X20002.2.17'] <- 'Self.Reported.Conditions.18'
names(ukb)[names(ukb) == 'X20002.2.18'] <- 'Self.Reported.Conditions.19'
names(ukb)[names(ukb) == 'X20002.2.19'] <- 'Self.Reported.Conditions.20'
names(ukb)[names(ukb) == 'X20002.2.20'] <- 'Self.Reported.Conditions.21'
names(ukb)[names(ukb) == 'X20002.2.21'] <- 'Self.Reported.Conditions.22'
names(ukb)[names(ukb) == 'X20002.2.22'] <- 'Self.Reported.Conditions.23'
names(ukb)[names(ukb) == 'X20002.2.23'] <- 'Self.Reported.Conditions.24'
names(ukb)[names(ukb) == 'X20002.2.24'] <- 'Self.Reported.Conditions.25'
names(ukb)[names(ukb) == 'X20002.2.25'] <- 'Self.Reported.Conditions.26'
names(ukb)[names(ukb) == 'X20002.2.26'] <- 'Self.Reported.Conditions.27'
names(ukb)[names(ukb) == 'X20002.2.27'] <- 'Self.Reported.Conditions.28'
names(ukb)[names(ukb) == 'X20002.2.28'] <- 'Self.Reported.Conditions.29'
names(ukb)[names(ukb) == 'X20002.2.29'] <- 'Self.Reported.Conditions.30'
names(ukb)[names(ukb) == 'X20002.2.30'] <- 'Self.Reported.Conditions.31'
names(ukb)[names(ukb) == 'X20002.2.31'] <- 'Self.Reported.Conditions.32'
names(ukb)[names(ukb) == 'X20002.2.32'] <- 'Self.Reported.Conditions.33'

# First step is to select only the samples that have imaging data
# 21003 column is 'Age when attended assesment center': 
#   .2.0 is the third visit which is the imaging visit.
#   If not NA then imaging was done. 
ukb_img <- ukb[!is.na(ukb$Age.Imaging.Visit),]
ukb_img <- ukb_img[!is.na(ukb_img$HV.Left),]
ukb_img <- ukb_img[!is.na(ukb_img$HV.Right),]

# Next we add link the imputed SNPs from the extracted raw files. 
# columns are: FID	IID	PAT	MAT	SEX	PHENOTYPE	[one or more SNP columns]
# we will select the IID and SNP columns from every file and join to full table by IID and eid
snp_2_raw <- read.csv("~/Desktop/ukb_data/filtered_chr2.raw" , header = TRUE , sep = "")
snp_5_raw <- read.csv("~/Desktop/ukb_data/filtered_chr5.raw" , header = TRUE , sep = "")
snp_7_raw <- read.csv("~/Desktop/ukb_data/filtered_chr7.raw" , header = TRUE , sep = "")
snp_9_raw <- read.csv("~/Desktop/ukb_data/filtered_chr9.raw" , header = TRUE , sep = "")
snp_12_raw <- read.csv("~/Desktop/ukb_data/filtered_chr12.raw" , header = TRUE , sep = "")

ukb_img <- merge(ukb_img , snp_2_raw[,c("IID","rs2268894_C")] , by.x="eid" , by.y="IID")
ukb_img <- merge(ukb_img , snp_5_raw[,c("IID","rs2289881_G")] , by.x="eid" , by.y="IID")
ukb_img <- merge(ukb_img , snp_7_raw[,c("IID","rs11979341_C")] , by.x="eid" , by.y="IID")
ukb_img <- merge(ukb_img , snp_9_raw[,c("IID","rs7020341_G")] , by.x="eid" , by.y="IID")
ukb_img <- merge(ukb_img , snp_12_raw[,c("IID","rs61921502_T","rs77956314_T")] , by.x="eid" , by.y="IID")

# calculate and add the PRS column
# RSID           Z-SCORE   FREQ     N     CHR  REF ALT   OUR_RSID
# rs77956314    -10.418   0.9160  26814   12    T   C    rs77956314_T
# rs61921502     9.017    0.8466  26814   12    T   G    rs61921502_T
# rs11979341    -6.755    0.6837  24484   7     C   G    rs11979341_C
# rs7020341      6.645    0.3590  26700   9     C   G    rs7020341_G     -- needs to be flipped
# rs2268894     -6.546    0.5412  26814   2     T   C    rs2268894_C     -- needs to be flipped
# rs2289881     -5.558    0.3544  26814   5     T   G    rs2289881_G     -- needs to be flipped

# Z-SCORES Are given in term of additive effects of REF
# So, to match the SNPs given in GWAS paper, we need to flip a few genotypes (the corresponding freqs)
# To flip the genotypes (0 -> 2 , 1 -> 1 , and 2 -> 0 ) perform 2 - Current.Genotype
ukb_img$rs7020341_C <- 2 - ukb_img$rs7020341_G
ukb_img$rs2268894_T <- 2 - ukb_img$rs2268894_C
ukb_img$rs2289881_T <- 2 - ukb_img$rs2289881_G

# Z-Scores are not the Beta effect sizes, though they can be a proxy
# We estimate the Beta from the Z-Score (z) , frequency (p), and sample size (n)
# using the formula: 
#   Beta = z / sqrt(2p(1− p)(n + z^2))
#   SE =1 / sqrt(2p(1− p)(n + z^2))
Z.Scores <- c(-10.418 , 9.017,  -6.755 ,6.645 , -6.546 , -5.558)
Freqs <- c( 0.9160 , 0.8466 , 0.6837 , (1-0.3590) , (1-0.5412) , (1-0.3544))
Ns <- c(26814 , 26814 , 24484 , 26700 , 26814 , 26814)
Betas <- Z.Scores / sqrt( 2*Freqs*(1-Freqs)*(Ns + Z.Scores^2) )
SEs <- 1 / sqrt( 2*Freqs*(1-Freqs)*(Ns + Z.Scores^2) )

ukb_img$PRS_6SNPS <- rowSums( ukb_img[,c("rs2268894_T","rs2289881_T","rs11979341_C","rs7020341_C","rs61921502_T","rs77956314_T")]) * sum(Betas)

ukb_img$PRS_6SNPS_ZSCORE <- rowSums( ukb_img[,c("rs2268894_T","rs2289881_T","rs11979341_C","rs7020341_C","rs61921502_T","rs77956314_T")]) * sum(Z.Scores)


# First PRS above showed no association so we used PRSice to find a PRS based on the full GWAS file and not just the 6 snps.
# Read in table that has the PRS caclculated at different p-value thresholds.
PRS_PRSICE <- read.table("~/Desktop/ukb_data/PRS-FINAL.all.score" , header = TRUE)
# rename columns before merging 
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X1e.08']  <- 'PRS_TH_1e.08'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X1e.07']  <- 'PRS_TH_1e.07'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X1e.06']  <- 'PRS_TH_1e.06'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X1e.05']  <- 'PRS_TH_1e.05'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.0001']  <- 'PRS_TH_1e.04'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.001']  <- 'PRS_TH_1e.03'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.01']  <- 'PRS_TH_0.01'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.05']  <- 'PRS_TH_0.05'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.1']  <- 'PRS_TH_0.1'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.2']  <- 'PRS_TH_0.2'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.4']  <- 'PRS_TH_0.4'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.5']  <- 'PRS_TH_0.5'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X0.75']  <- 'PRS_TH_0.75'
names(PRS_PRSICE)[names(PRS_PRSICE) == 'X1']  <- 'PRS_TH_1'

ukb_img <- merge(ukb_img , PRS_PRSICE , by.x="eid" , by.y = "IID" , all.x=TRUE )
# Exclusions

# Need to exclude a few sample EID's
private_eids <- c('1090030','1091961','1230112','1325363','1503907','1546281','1779627','2072543','2129742','2328077','2381848',
                  '2562577','2682290','2850777','2991997','3279656','3303963','3629688','3640543','3819653','3850604','3874773',
                  '3983605','3999481','4217643','4278871','4362486','4588727','4652426','4684825','4838503','5122856','5203595',
                  '5494205','5708826','5912074','5954147')

ukb_img_ex <- ukb_img[! ukb_img$eid %in% private_eids , ]

# Need to include only White/British ethnicity.
# Column 21000 is "Ethnic Background" and is the self-reported ethnicity of sample.
# value 1001 means "British". So we are only including those for now.
ukb_img_ex <- ukb_img_ex[ukb_img_ex$Ethnic.Background==1001,]

# Need to exclude people with Neurological or Psychiatric Disorders
# and people with substance abuse or history of head trauma
# and people with cardiovascular disorders
# 1079  cardiomyopathy
# 1240  neurological injury/trauma
# 1243  psychological/psychiatric problem
# 1258  chronic/degenerative neurological problem
# 1261	multiple sclerosis
# 1262	parkinsons disease
# 1263	dementia/alzheimers/cognitive impairment
# 1264	epilepsy
# 1266	head injury
# 1289	schizophrenia
# 1408  alcohol dependency
# 1409  opioid dependency
# 1410  other substance abuse/dependency
# 1425	cerebral aneurysm
# 1434  other neurological problem
# 1491	brain haemorrhage
exclude_conditions <- c( 1079 , 1240 , 1243 , 1258 , 1261 , 1262 , 1263 , 1264 , 1266 , 1289 , 
                         1408 , 1409 , 1410 , 1425 , 1434 , 1491 )
excluded_rows <- ukb_img_ex
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.1  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.2  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.3  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.4  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.5  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.6  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.7  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.8  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.9  %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.10 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.11 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.12 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.13 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.14 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.15 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.16 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.17 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.18 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.19 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.20 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.21 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.22 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.23 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.24 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.25 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.26 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.27 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.28 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.29 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.30 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.31 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.32 %in% exclude_conditions , ]
ukb_img_ex <- ukb_img_ex[! ukb_img_ex$Self.Reported.Conditions.33 %in% exclude_conditions , ]

excluded_rows <- subset( excluded_rows , !( excluded_rows$eid %in% ukb_img_ex$eid) )

# Second step is to filter out outliers.
# This is done by filtering out all volumes that are more than 5 MAE's away from the mean. 
# 25019.2.0 is left HV, 25020.2.0 is right HV.
# Find the mean of both HV
ukb_img_ex_hv_left_mean  <- mean(ukb_img_ex$HV.Left , na.rm = TRUE)
ukb_img_ex_hv_right_mean <- mean(ukb_img_ex$HV.Right , na.rm = TRUE)
# Then the MAE: 1/n * sum( x - mean(x) )
ukb_img_ex_hv_left_mae  <- sum(abs(ukb_img_ex$HV.Left - ukb_img_ex_hv_left_mean ) , na.rm = TRUE)  / sum(!is.na(ukb_img_ex$HV.Left))
ukb_img_ex_hv_right_mae <- sum(abs(ukb_img_ex$HV.Right - ukb_img_ex_hv_right_mean ) , na.rm = TRUE) / sum(!is.na(ukb_img_ex$HV.Right))
# Do the outlier filtering
ukb_img_ex_outliers <- ukb_img_ex[ abs(ukb_img_ex$HV.Left - ukb_img_ex_hv_left_mean) < 5 * ukb_img_ex_hv_left_mae , ]
ukb_img_ex_outliers <- ukb_img_ex_outliers[ abs(ukb_img_ex$HV.Right - ukb_img_ex_hv_right_mean) < 5 * ukb_img_ex_hv_right_mae , ]



# Third step is to correct for some variables by regressing them out. 
# Head size indicator which in this case is:
# 25000.0.0 'Volumetric scaling from T1 head image to standard space'
# Scan Date: 
# 53.2.0 'Scan Date'

mod_left  <- lm( HV.Left  ~  Volumetric.scaling.from.T1.head.image.to.standard.space 
                           + as.numeric(as.Date(Scan.Date)) ,
                           # + Sex , # + Age.Imaging.Visit
                           ukb_img_ex_outliers , na.action = na.exclude )

mod_right <- lm( HV.Right ~ Volumetric.scaling.from.T1.head.image.to.standard.space 
                          + as.numeric(as.Date(Scan.Date)) ,
                          # + Sex , # + Age.Imaging.Visit
                          ukb_img_ex_outliers , na.action = na.exclude )

mean_scaling_factor <- mean(ukb_img_ex_outliers$Volumetric.scaling.from.T1.head.image.to.standard.space , na.rm = TRUE)
mean_scan_date <- mean(as.numeric(as.Date(ukb_img_ex_outliers$Scan.Date)) , na.rm = TRUE)

# Now shift the residuals (centered around zero) so thet their values are back in the normal hv volume range. 
# we do this with the equation:  Corrected_volume = Intercept + Residuals + ( Beta of ICV * mean ICV)
clean_hv_volumes_left <- coefficients(mod_left)[1] + resid(mod_left) + (coefficients(mod_left)[2] * mean_scaling_factor) + (coefficients(mod_left)[3] * mean_scan_date)
clean_hv_volumes_right <- coefficients(mod_right)[1] + resid(mod_right) + (coefficients(mod_right)[2] * mean_scaling_factor) + (coefficients(mod_right)[3] * mean_scan_date)

ukb_img_ex_outliers$clean_hv_volumes_left <- clean_hv_volumes_left
ukb_img_ex_outliers$clean_hv_volumes_right <- clean_hv_volumes_right


##### PERFORM WINDOW ANALYSIS #########

# Now we will generate a table with binned information to begin constructing the nomogram.
# The bins should each contain 10% of the samples

# sort by age
ukb_img_ex_outliers <- ukb_img_ex_outliers[order(ukb_img_ex_outliers$Age.Imaging.Visit),]

# stratify by gender
ukb_img_ex_outliers_males   <- ukb_img_ex_outliers[ukb_img_ex_outliers$Sex == 1,]
ukb_img_ex_outliers_females <- ukb_img_ex_outliers[ukb_img_ex_outliers$Sex == 0,]

# call the function that does the window analysis
bins_male <- WINDOW_ANALYSIS( ukb_img_ex_outliers_males )
bins_female <- WINDOW_ANALYSIS( ukb_img_ex_outliers_females )

##### PLOTTING  #######

PLOT_NOMOGRAM( bins_male , 20 , 1 , "Nomogram:\n Males Left HV")
PLOT_NOMOGRAM( bins_male , 20 , 0 , "Nomogram:\n Males Right HV")
PLOT_NOMOGRAM( bins_female , 20 , 1 , "Nomogram:\n Females Left HV")
PLOT_NOMOGRAM( bins_female , 20 , 0 , "Nomogram:\n Females Right HV")


###### GENETIC ANALYSIS #######

par(mfrow=c(3,3))
for( GENDER in c("Male","Female","Both Genders")){
  if (GENDER=="Male") table <- ukb_img_ex_outliers_males
  if (GENDER=="Female") table <- ukb_img_ex_outliers_females
  if (GENDER=="Both Genders") table <- ukb_img_ex_outliers
  
  for(HEMISPHERE in c("Left","Right","Both Hemispheres")){
    if (HEMISPHERE=="Left") hv <- table$clean_hv_volumes_left
    if (HEMISPHERE=="Right") hv <- table$clean_hv_volumes_right
    if (HEMISPHERE=="Both Hemispheres") hv <- table$clean_hv_volumes_left + table$clean_hv_volumes_right
    
    HV_SNP <- matrix(c(0,0,0,0),ncol=4,nrow=14)
    colnames(HV_SNP) <- c("Slope" , "Range" , "P-Value" , "R-Squared")
    
    thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.06','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                    'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')
    rownames(HV_SNP) <- thresholds
    i=1
    for( col_name in thresholds ){
      x <- lm( hv ~ table[,col_name] )
      HV_SNP[i,1] <- summary(x)$coefficients[2,1]
      HV_SNP[i,2] <- max(table[,col_name] , na.rm=TRUE) - min(table[,col_name] , na.rm=TRUE)
      HV_SNP[i,3] <- summary(x)$coefficients[2,4]
      HV_SNP[i,4] <- summary(x)$r.squared
      i<-i+1
    }
    
    bplot <- barplot(main = paste(GENDER , " " , HEMISPHERE) , 
                     HV_SNP[,4] , xlab = "P-Value Threshold" , ylab = "R-Squared" , ylim = c(0,0.014) , 
            names=c("1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.05","0.1" , "0.2" , "0.4" , "0.5" , "0.75" , "1"))
    
    text( x=c(1:14)*1.2,  y=HV_SNP[,4]+0.0001, formatC(HV_SNP[,3], format = "e", digits = 2) , srt=80 , adj=c(0,-1))
  }
}

