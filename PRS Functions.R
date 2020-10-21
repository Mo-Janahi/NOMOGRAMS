# File Containing all functions that deal with Polygenic Risk Score
# By: Mohammed Janahi
#
#


MAKE_PRS_PLOTS <- function( table_male , table_female , confounder_columns , 
                            hv_l_column = "clean_hv_left", hv_r_column = "clean_hv_right",
                            hv_bl_column = "clean_hv_bilateral"){
  
  full <- merge(table_male , table_female , by=names(table_male) , all.x = TRUE , all.y = TRUE)
  # show what the prs percentiles look like 
  par(mfrow=c(3,6))
  for( GENDER in c("Male","Female","Both Genders")){
    if (GENDER=="Male")         table <- table_male
    if (GENDER=="Female")       table <- table_female
    if (GENDER=="Both Genders") table <- full
    for(HEMISPHERE in c("Left","Right","Both Hemispheres")){
      if (HEMISPHERE=="Left")             hv <- hv_l_column
      if (HEMISPHERE=="Right")            hv <- hv_r_column
      if (HEMISPHERE=="Both Hemispheres") hv <- hv_bl_column
      
      HV_SNP <- matrix(c(0,0,0,0),ncol=4,nrow=14)
      colnames(HV_SNP) <- c("Slope" , "Range" , "P-Value" , "R-Squared")
      
      thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.06','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                      'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')
      rownames(HV_SNP) <- thresholds
      i=1
      for( col_name in thresholds ){
        x <- lm( paste( hv , " ~ ",col_name," + " , paste( confounder_columns , collapse = " + ") ) , data = table )
        HV_SNP[i,1] <- summary(x)$coefficients[2,1]
        HV_SNP[i,2] <- max(table[,col_name] , na.rm=TRUE) - min(table[,col_name] , na.rm=TRUE)
        HV_SNP[i,3] <- summary(x)$coefficients[2,4]
        HV_SNP[i,4] <- summary(x)$r.squared
        i<-i+1
      }
      
      bplot <- barplot(main = paste(GENDER , " " , HEMISPHERE) , xpd = FALSE , 
                       HV_SNP[,4] , xlab = "P-Value Threshold" , ylab = "R-Squared" , ylim = c(summary(HV_SNP[,4])[1]-0.001,summary(HV_SNP[,4])[6]+0.01) , 
                       names=c("1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","0.01","0.05","0.1" , "0.2" , "0.4" , "0.5" , "0.75" , "1"))
      
      text( x=c(1:14)*1.2,  y=HV_SNP[,4]+0.0001, formatC(HV_SNP[,3], format = "e", digits = 2) , srt=80 , adj=c(0,-1))
      
      max_col_name <-  names(HV_SNP[,4])[match(max(HV_SNP[,4]),HV_SNP[,4])]
      hv <- table[ , hv]
      hv <- hv[  order(table[,max_col_name]) ]
      
      means <- 1:100
      for ( i in 1:100 ){
        # THIS METHOD CAlCULATESS MEANS AS PERCENTILE OF SAMPLES
        start_ind = nrow(table) * ((i-1)/100)
        end_ind = nrow(table) * (i/100)
        window <- hv[start_ind:end_ind]
        means[i] <- mean((window))
      }
      plot(means , main = paste("Percentiles for Threshold: " , max_col_name) )
      
    }
  }
  return(HV_SNP)
}

