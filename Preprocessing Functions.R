



REGRESS_OUT_COLUMN <- function( col_name , regress_column , regress_column_2 = NA , regress_column_3 = NA){
  
  if ( is.na(regress_column_2[1]) ) regress_column_2 <- regress_column * 0
  if ( is.na(regress_column_3[1]) ) regress_column_3 <- regress_column * 0
  
  mod  <- lm( col_name  ~  regress_column + regress_column_2 + regress_column_3, na.action = na.exclude )
  
  mean_factor <- mean( regress_column , na.rm = TRUE)
  
  # Now shift the residuals (centered around zero) so thet their values are back in the normal hv volume range. 
  # we do this with the equation:  Corrected_volume = Intercept + Residuals + ( Beta of ICV * mean ICV)
  clean_col <- coefficients(mod)[1] + resid(mod) + (coefficients(mod)[2] * mean_factor)
  
  return(clean_col)
}

REGRESS_OUT <- function( hv_col_name , ICV_col_name , table ){
  
  # regressing variables out. 
  mod  <- lm( table[,hv_col_name]  ~  table[,ICV_col_name]
              + as.numeric(as.Date(table[,"Scan.Date"])) ,
              na.action = na.exclude )
  
  mean_scaling_factor <- mean(table[,ICV_col_name] , na.rm = TRUE)
  mean_scan_date <- mean(as.numeric(as.Date(table[,"Scan.Date"])) , na.rm = TRUE)
  
  # Now shift the residuals (centered around zero) so thet their values are back in the normal hv volume range. 
  # we do this with the equation:  Corrected_volume = Intercept + Residuals + ( Beta of ICV * mean ICV)
  clean_hv <- coefficients(mod)[1] + resid(mod) + (coefficients(mod)[2] * mean_scaling_factor) + (coefficients(mod)[3] * mean_scan_date)
  
  return(clean_hv)
}

READ_PRS_TABLE <- function( table_path , prefix=NA){
  # Read in table that has the PRS calculated at different p-value thresholds.
  #PRS <- read.table("~/Desktop/ukb_data/PRS-FINAL.all.score" , header = TRUE)
  #PRS <- read.table("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/PRS_CLN_HV_40k.all.score" , header = TRUE)
  PRS <- read.table(table_path , header = TRUE)
  
  # rename columns before merging 
  names(PRS)[names(PRS) == 'X1e.08']  <- ifelse( is.na(prefix) , 'PRS_TH_1e.08' , paste(prefix,'PRS_TH_1e.08',sep = "_") )
  names(PRS)[names(PRS) == 'X1e.07']  <- ifelse( is.na(prefix) , 'PRS_TH_1e.07' , paste(prefix,'PRS_TH_1e.07',sep = "_") )
  names(PRS)[names(PRS) == 'X1e.06']  <- ifelse( is.na(prefix) , 'PRS_TH_1e.06' , paste(prefix,'PRS_TH_1e.06',sep = "_") )
  names(PRS)[names(PRS) == 'X1e.05']  <- ifelse( is.na(prefix) , 'PRS_TH_1e.05' , paste(prefix,'PRS_TH_1e.05',sep = "_") )
  names(PRS)[names(PRS) == 'X0.0001'] <- ifelse( is.na(prefix) , 'PRS_TH_1e.04' , paste(prefix,'PRS_TH_1e.04',sep = "_") )
  names(PRS)[names(PRS) == 'X0.001']  <- ifelse( is.na(prefix) , 'PRS_TH_1e.03' , paste(prefix,'PRS_TH_1e.03',sep = "_") )
  names(PRS)[names(PRS) == 'X0.01']   <- ifelse( is.na(prefix) , 'PRS_TH_0.01'  , paste(prefix,'PRS_TH_0.01',sep = "_") )
  names(PRS)[names(PRS) == 'X0.05']   <- ifelse( is.na(prefix) , 'PRS_TH_0.05'  , paste(prefix,'PRS_TH_0.05',sep = "_") )
  names(PRS)[names(PRS) == 'X0.1']    <- ifelse( is.na(prefix) , 'PRS_TH_0.1'   , paste(prefix,'PRS_TH_0.1',sep = "_") )
  names(PRS)[names(PRS) == 'X0.2']    <- ifelse( is.na(prefix) , 'PRS_TH_0.2'   , paste(prefix,'PRS_TH_0.2',sep = "_") )
  names(PRS)[names(PRS) == 'X0.4']    <- ifelse( is.na(prefix) , 'PRS_TH_0.4'   , paste(prefix,'PRS_TH_0.4',sep = "_") )
  names(PRS)[names(PRS) == 'X0.5']    <- ifelse( is.na(prefix) , 'PRS_TH_0.5'   , paste(prefix,'PRS_TH_0.5',sep = "_") )
  names(PRS)[names(PRS) == 'X0.75']   <- ifelse( is.na(prefix) , 'PRS_TH_0.75'  , paste(prefix,'PRS_TH_0.75',sep = "_") )
  names(PRS)[names(PRS) == 'X1']      <- ifelse( is.na(prefix) , 'PRS_TH_1'     , paste(prefix,'PRS_TH_1',sep = "_") )
  
  return(PRS)
  
  # for reference 
  ## OUR FIRST HYPOTHESIS  ##
  # First thing we tested was whether 6 snps identified in another paper for HV
  # Next we add link the imputed SNPs from the extracted raw files. 
  # columns are: FID	IID	PAT	MAT	SEX	PHENOTYPE	[one or more SNP columns]
  # we will select the IID and SNP columns from every file and join to full table by IID and eid
  #snp_2_raw <- read.csv("~/Desktop/ukb_data/filtered_chr2.raw" , header = TRUE , sep = "")
  #snp_5_raw <- read.csv("~/Desktop/ukb_data/filtered_chr5.raw" , header = TRUE , sep = "")
  #snp_7_raw <- read.csv("~/Desktop/ukb_data/filtered_chr7.raw" , header = TRUE , sep = "")
  #snp_9_raw <- read.csv("~/Desktop/ukb_data/filtered_chr9.raw" , header = TRUE , sep = "")
  #snp_12_raw <- read.csv("~/Desktop/ukb_data/filtered_chr12.raw" , header = TRUE , sep = "")
  
  #ukb_img <- merge(ukb_img , snp_2_raw[,c("IID","rs2268894_C")] , by.x="eid" , by.y="IID")
  #ukb_img <- merge(ukb_img , snp_5_raw[,c("IID","rs2289881_G")] , by.x="eid" , by.y="IID")
  #ukb_img <- merge(ukb_img , snp_7_raw[,c("IID","rs11979341_C")] , by.x="eid" , by.y="IID")
  #ukb_img <- merge(ukb_img , snp_9_raw[,c("IID","rs7020341_G")] , by.x="eid" , by.y="IID")
  #ukb_img <- merge(ukb_img , snp_12_raw[,c("IID","rs61921502_T","rs77956314_T")] , by.x="eid" , by.y="IID")
  
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
  #ukb_img$rs7020341_C <- 2 - ukb_img$rs7020341_G
  #ukb_img$rs2268894_T <- 2 - ukb_img$rs2268894_C
  #ukb_img$rs2289881_T <- 2 - ukb_img$rs2289881_G
  
  # Z-Scores are not the Beta effect sizes, though they can be a proxy
  # We estimate the Beta from the Z-Score (z) , frequency (p), and sample size (n)
  # using the formula: 
  #   Beta = z / sqrt(2p(1− p)(n + z^2))
  #   SE =1 / sqrt(2p(1− p)(n + z^2))
  #Z.Scores <- c(-10.418 , 9.017,  -6.755 ,6.645 , -6.546 , -5.558)
  #Freqs <- c( 0.9160 , 0.8466 , 0.6837 , (1-0.3590) , (1-0.5412) , (1-0.3544))
  #Ns <- c(26814 , 26814 , 24484 , 26700 , 26814 , 26814)
  #Betas <- Z.Scores / sqrt( 2*Freqs*(1-Freqs)*(Ns + Z.Scores^2) )
  #SEs <- 1 / sqrt( 2*Freqs*(1-Freqs)*(Ns + Z.Scores^2) )
  
  #ukb_img$PRS_6SNPS <- rowSums( ukb_img[,c("rs2268894_T","rs2289881_T","rs11979341_C","rs7020341_C","rs61921502_T","rs77956314_T")]) * sum(Betas)
  
  #ukb_img$PRS_6SNPS_ZSCORE <- rowSums( ukb_img[,c("rs2268894_T","rs2289881_T","rs11979341_C","rs7020341_C","rs61921502_T","rs77956314_T")]) * sum(Z.Scores)
  
  # First PRS above showed no association with HV 
  # So we used PRSice to find a PRS based on the full GWAS file and not just the 6 snps.
}

MAD_FILTER <- function( table , columns , threshold){
  
  for (col in columns) {
    vals <- table[,col]
    mean<- mean(vals , na.rm = TRUE)
    mae <- sum(abs(vals - mean ) , na.rm = TRUE)  / sum(!is.na(vals))
    table <- table[ which( abs(vals- mean) < 5 * mae) , ]
  }
  return (table)
  
  #ukb_img_ex_hv_left_mean  <- mean(ukb_img_ex$HV.Left , na.rm = TRUE)
  # Then the MAE: 1/n * sum( x - mean(x) )
  #ukb_img_ex_hv_left_mae  <- sum(abs(ukb_img_ex$HV.Left - ukb_img_ex_hv_left_mean ) , na.rm = TRUE)  / sum(!is.na(ukb_img_ex$HV.Left))
  
  # Repeat for right hemisphere
  #ukb_img_ex_hv_right_mean <- mean(ukb_img_ex$HV.Right , na.rm = TRUE)
  #ukb_img_ex_hv_right_mae <- sum(abs(ukb_img_ex$HV.Right - ukb_img_ex_hv_right_mean ) , na.rm = TRUE) / sum(!is.na(ukb_img_ex$HV.Right))
  
  # Do the outlier filtering
  #ukb_img_ex_outliers <- ukb_img_ex[ which( (abs(ukb_img_ex$HV.Right - ukb_img_ex_hv_right_mean) < 5 * ukb_img_ex_hv_right_mae)
  #                                          | (abs(ukb_img_ex$HV.Left - ukb_img_ex_hv_left_mean) < 5 * ukb_img_ex_hv_left_mae) ) , ]
  
  
}

PREPROCESS_UKBB <- function( ukb ){
  
  #ukb <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200529.csv")
  #ukb_ext <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200602_FSA.csv")
  # Merge the two tables
  #ukb <- merge(ukb , ukb_ext , by=intersect(names(ukb),names(ukb_ext) )) # merge by intersect makes sure we don't duplicate columns
  
  # Andre put together a file with more accurate ages. Ages in table are just years without months/decimals.
  # but taking exact birth date and scan date into account will give more accurate ages. 
  ukb_ages <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/ukb_all_subjs.csv")
  # create a column with the age from the latest imaging visit (some subjects came for multiple visits)
  ukb_ages$AGE_Latest <- ukb_ages$AGE2
  ukb_ages$AGE_Latest [is.na(ukb_ages$AGE2)] <- ukb_ages$AGE[is.na(ukb_ages$AGE2)]
  
  ukb <- merge(ukb , ukb_ages , by="eid") 
  
  
  # For readability, rename some columns
  names(ukb)[names(ukb) == 'X31.0.0']    <- 'Sex'
  names(ukb)[names(ukb) == 'X21003.0.0'] <- 'Age.When.Attended.Assesment.Center'
  names(ukb)[names(ukb) == 'X21003.2.0'] <- 'Age.Imaging.Visit'
  names(ukb)[names(ukb) == 'X25000.2.0'] <- 'Volumetric.scaling.from.T1.head.image.to.standard.space'
  names(ukb)[names(ukb) == 'X53.2.0']    <- 'Scan.Date'
  names(ukb)[names(ukb) == 'X25019.2.0'] <- 'HV.Left.FSL'
  names(ukb)[names(ukb) == 'X25020.2.0'] <- 'HV.Right.FSL'
  names(ukb)[names(ukb) == 'X26562.2.0'] <- 'HV.Left'
  names(ukb)[names(ukb) == 'X26593.2.0'] <- 'HV.Right'
  names(ukb)[names(ukb) == 'X21000.0.0'] <- 'Ethnic.Background'
  for( i in 1:length(grep("X22009.0." , names(ukb))) )
    names(ukb)[names(ukb) == paste('X22009.0.',i,sep="")] <- paste('Genetic.PC.',i,sep="")
  for( i in 0:length(grep("X20002.2." , names(ukb))) )
    names(ukb)[names(ukb) == paste('X20002.2.',i,sep="")]  <- paste('Self.Reported.Conditions.',(i+1),sep="")
  
  # First step is to select only the samples that have imaging data
  ukb_img <- ukb[!is.na(ukb$Age.Imaging.Visit),]
  ukb_img <- ukb_img[!is.na(ukb_img$HV.Left),]
  ukb_img <- ukb_img[!is.na(ukb_img$HV.Right),]
  
  
  # Exclusions
  
  # Need to exclude a few sample EID's (patient choose to be excluded)
  private_eids <- c('1090030','1091961','1230112','1325363','1503907','1546281','1779627','2072543','2129742','2328077','2381848',
                    '2562577','2682290','2850777','2991997','3279656','3303963','3629688','3640543','3819653','3850604','3874773',
                    '3983605','3999481','4217643','4278871','4362486','4588727','4652426','4684825','4838503','5122856','5203595',
                    '5494205','5708826','5912074','5954147')
  
  ukb_img_ex <- ukb_img[! ukb_img$eid %in% private_eids , ]
  
  # Need to include only White/British ethnicity.
  # value 1001 means "British". So we are only including those for now.
  ukb_img_ex <- ukb_img_ex[which(ukb_img_ex$Ethnic.Background==1001),]
  
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
  
  for( i in 1:34 )
    ukb_img_ex <- ukb_img_ex[! ukb_img_ex[,paste("Self.Reported.Conditions.",i,sep="")] %in% exclude_conditions , ]
  
  excluded_rows <- subset( excluded_rows , !( excluded_rows$eid %in% ukb_img_ex$eid) )
  
  # stratify by gender
  ukb_img_ex_outliers_male   <- ukb_img_ex[ which(ukb_img_ex$Sex == 1),]
  ukb_img_ex_outliers_female <- ukb_img_ex[ which(ukb_img_ex$Sex == 0),]
  
  
  # Outliers
  
  # next step is to filter out outliers.
  # This is done by filtering out all volumes that are more than 5 MAE's away from the mean. 
  # Find the mean of both HV
  ukb_img_ex_outliers_male <- MAD_FILTER(ukb_img_ex_outliers_male , c("HV.Left","HV.Right") , 5)
  ukb_img_ex_outliers_female <- MAD_FILTER(ukb_img_ex_outliers_female , c("HV.Left","HV.Right") , 5)
  
  
  
  # Confounders
  
  # next step is to correct for some variables by regressing them out. 
  
  # do the regressing out of variables independently for each gender 
  ukb_img_ex_outliers_male$clean_hv_left <- REGRESS_OUT("HV.Left" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_male)
  ukb_img_ex_outliers_male$clean_hv_right <- REGRESS_OUT("HV.Right" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_male)
  ukb_img_ex_outliers_male$clean_hv_bilateral <- ( ukb_img_ex_outliers_male$clean_hv_left + ukb_img_ex_outliers_male$clean_hv_right )/2
  
  ukb_img_ex_outliers_female$clean_hv_left <- REGRESS_OUT("HV.Left" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_female)
  ukb_img_ex_outliers_female$clean_hv_right <- REGRESS_OUT("HV.Right" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_female)
  ukb_img_ex_outliers_female$clean_hv_bilateral <- (ukb_img_ex_outliers_female$clean_hv_left + ukb_img_ex_outliers_female$clean_hv_right)/2

  
  # Genetics
  
  
  # If we want PRSICE to also generate PRS bar plots, then
  # we need to export the HV in a file with 3 columns:
  # FID : the eid from our table
  # IID : again the eid from our table
  # HV : corrected bilateral hippocampal volume
  
  #full <- merge(ukb_img_ex_outliers_male , ukb_img_ex_outliers_female , by=names(ukb_img_ex_outliers_male) , all.x = TRUE , all.y = TRUE)
  #to_save <- full[ , c("eid" , "eid" , "clean_hv_bilateral")]
  #names(to_save) <- c("FID" , "IID" , "clean_hv_bilateral")
  #write.table( to_save , "/Users/mohammedjanahi/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_cal_merged_maf01.pheno", 
  #             append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE , quote = FALSE)
  
  # this table goes to PRSice where a PRS is calculated with the command:
  
  # -------->
  #         Rscript ../PRSice_mac/PRSice.R --prsice ./../PRSice_mac/PRSice \
  #                    --base ../CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL \
  #                    --target ../ukb_cal_merged_maf01 --missing SET_ZERO \
  #                    --geno 0.1 --maf 0.05 --no-regress --fastscore --all-score --score sum \
  #                    --binary-target F --thread 1 \
  #                    --beta --stat Beta --A1 Allele1 --snp RSNUMBERS --pvalue P.value \
  #                    --bar-levels 0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.2,0.4,0.5,0.75,1.0 \
  #                    --out PRS_HV_FREESURFER_40k
  #
  # <--------
  # then we read in the results of this call for analysis
  
  # PRS_PRSICE <- READ_PRS_TABLE("~/Desktop/ukb_data/PRS-FINAL.all.score")
  # PRS_PRSICE <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/PRS_CLN_HV_40k.all.score")
  HV_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/PRS_HV_FREESURFER_40k.all.score")
  ukb_img_ex_outliers_male <- merge(ukb_img_ex_outliers_male , HV_PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )
  ukb_img_ex_outliers_female <- merge(ukb_img_ex_outliers_female , HV_PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )
  
  AD_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/PRSice_mac/AD_PRS_HV_40k.all.score" , "AD")
  ukb_img_ex_outliers_male <- merge(ukb_img_ex_outliers_male , AD_PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )
  ukb_img_ex_outliers_female <- merge(ukb_img_ex_outliers_female , AD_PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )

  ICV_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/PRSice_mac/ICV_PRS_HV_40k.all.score" , "ICV")
  ukb_img_ex_outliers_male <- merge(ukb_img_ex_outliers_male , ICV_PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )
  ukb_img_ex_outliers_female <- merge(ukb_img_ex_outliers_female , ICV_PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )
  
  SE_0.01 <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/PRSice_mac/PRS_HV_SE_0.01.all.score" , "SE_0.01")
  ukb_img_ex_outliers_male <- merge(ukb_img_ex_outliers_male , SE_0.01 , by.x="eid" , by.y = "IID" , all.x=TRUE )
  ukb_img_ex_outliers_female <- merge(ukb_img_ex_outliers_female , SE_0.01 , by.x="eid" , by.y = "IID" , all.x=TRUE )
  
  #
  SE_0.0087 <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/PRSice_mac/PRS_HV_SE_0.0087.all.score" , "SE_0.0087")
  ukb_img_ex_outliers_male <- merge(ukb_img_ex_outliers_male , SE_0.0087 , by.x="eid" , by.y = "IID" , all.x=TRUE )
  ukb_img_ex_outliers_female <- merge(ukb_img_ex_outliers_female , SE_0.0087 , by.x="eid" , by.y = "IID" , all.x=TRUE )
  
  
  return( list(ukb_img_ex_outliers_male = ukb_img_ex_outliers_male,
               ukb_img_ex_outliers_female = ukb_img_ex_outliers_female) )
}

PREPROCESS_ADNI <- function( ADNI_table ){
   
  # rename some columns for readability
  # col names : 
  # ST88SV : Volume (WM Parcellation) of RightHippocampus
  # ST29SV : Volume (WM Parcellation) of LeftHippocampus
  names(ADNI_table)[names(ADNI_table) == 'ST88SV'] <- 'HV_Right'
  names(ADNI_table)[names(ADNI_table) == 'ST29SV'] <- 'HV_Left'
  ADNI_table$Sex <- (ADNI_table$PTGENDER == "Male")
  
  
  ADNI_Gentic_PCs <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/adni_genetic_pcs.csv" , header = TRUE)
  
  ADNI_table <- merge(ADNI_table , ADNI_Gentic_PCs , by.x = "RID.x" , by.y = "RID" )
  
  #Rscript ../ukb_data/PRSice_mac/PRSice.R --prsice ./../ukb_data/PRSice_mac/PRSice \
  #          --base ../ukb_data/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL \
  #          --target merged_imputed_maf001_geno010_4batches \
  #           --extract ../ukb_data/with_imaging_40k/ukb_cal_merged_maf01.rsids \
  #          --missing SET_ZERO \
  #          --geno 0.1 --maf 0.05 --no-regress --fastscore --all-score --score sum \
  #          --binary-target F --thread 1 \
  #          --beta --stat Beta --A1 Allele1 --snp RSNUMBERS --pvalue P.value \
  #          --bar-levels 0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.2,0.4,0.5,0.75,1.0 \
  #          --out PRS_HV_FREESURFER_ADNI
  
  #ADNI_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/PRS_HV_FREESURFER_ADNI.all.score")
  #ADNI_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/PRS_HV_FREESURFER_ADNI3.all.score")
  ADNI_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/PRS_HV_FREESURFER_ADNI_NO_EXTRACT.all.score")
  
  ADNI_table <- merge(ADNI_table , ADNI_PRS , by.x="PTID" , by.y="IID" , all.x=TRUE)
  #ADNI_table <- merge(ADNI_table , ADNI_PRS , by.x="PTID" , by.y="IID")
  
  # We need to clean up the data before we can use it
  # exclude non-british/white samples 
  ADNI_filter <- ADNI_table[ ADNI_table$PTETHCAT == "Not Hisp/Latino" , ]
  ADNI_filter <- ADNI_filter[ ADNI_filter$PTRACCAT == "White" , ]
  
  ADNI_filter$TRUE.AGE <- ADNI_filter$AGE + ADNI_filter$Years.bl
  
  ADNI_filter[ ADNI_filter$RHIPQC == "Fail" , "HV_Right" ] = NA
  ADNI_filter[ ADNI_filter$LHIPQC == "Fail" , "HV_Left" ] = NA
  
  # Male strata
  ADNI_filter_male <- ADNI_filter[ ADNI_filter$PTGENDER == "Male", ]
  
  # correct for ICV and scan date
  ADNI_filter_male$clean_hv_right <- REGRESS_OUT_COLUMN(ADNI_filter_male$HV_Right , ADNI_filter_male$ICV) #,
#                                                        as.numeric(ADNI_filter_male$EXAMDATE.x))
  ADNI_filter_male$clean_hv_left <- REGRESS_OUT_COLUMN(ADNI_filter_male$HV_Left , ADNI_filter_male$ICV) #,
 #                                                      as.numeric(ADNI_filter_male$EXAMDATE.x))
  
  ADNI_filter_male$clean_hv_bilateral <- ( ADNI_filter_male$clean_hv_left + ADNI_filter_male$clean_hv_right ) / 2
  
  # Female strata
  ADNI_filter_female <- ADNI_filter[ ADNI_filter$PTGENDER == "Female", ]
  # correct for ICV
  ADNI_filter_female$clean_hv_right <- REGRESS_OUT_COLUMN(ADNI_filter_female$HV_Right , ADNI_filter_female$ICV )#,
#                                                          as.numeric(ADNI_filter_female$EXAMDATE.x))
  ADNI_filter_female$clean_hv_left <- REGRESS_OUT_COLUMN(ADNI_filter_female$HV_Left , ADNI_filter_female$ICV )#,
#                                                         as.numeric(ADNI_filter_female$EXAMDATE.x))
  ADNI_filter_female$clean_hv_bilateral <- ( ADNI_filter_female$clean_hv_left + ADNI_filter_female$clean_hv_right ) / 2
  
  return( list( ADNI_male = ADNI_filter_male ,
                ADNI_female = ADNI_filter_female ) )
}