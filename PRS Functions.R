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
  #par(mfrow=c(1,2))
  for( GENDER in c("Male","Female","Both Sexes")){
  #for( GENDER in c("Both Sexes")){
    if (GENDER=="Male")         table <- table_male
    if (GENDER=="Female")       table <- table_female
    if (GENDER=="Both Sexes") table <- full
    for(HEMISPHERE in c("Left","Right","Both Hemispheres")){
    #for(HEMISPHERE in c("Both Hemispheres")){
      if (HEMISPHERE=="Left")             hv <- hv_l_column
      if (HEMISPHERE=="Right")            hv <- hv_r_column
      if (HEMISPHERE=="Both Hemispheres") hv <- hv_bl_column
      
      #table[,hv] <- (table[,hv] - min(table[,hv],na.rm = TRUE)) / (max(table[,hv],na.rm = TRUE) - min(table[,hv],na.rm = TRUE))
      
      thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.06','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                      'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')
      
      thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                      'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')
      
      #thresholds <- c('PRS_TH_0.75','PRS_TH_0.75')
      
      #thresholds <- paste("ICV",thresholds,sep="_")
      
      HV_SNP <- matrix(c(0,0,0,0),ncol=4,nrow=length(thresholds))
      
      colnames(HV_SNP) <- c("Slope" , "Range" , "P-Value" , "R-Squared")
      rownames(HV_SNP) <- thresholds
      i=1
      for( col_name in thresholds ){
        #table[,col_name] <- (table[,col_name] - min(table[,col_name],na.rm = TRUE)) / (max(table[,col_name],na.rm = TRUE) - min(table[,col_name],na.rm = TRUE))
        #table[,col_name] <- (table[,hv] - min(table[,hv],na.rm = TRUE)) / (max(table[,hv],na.rm = TRUE) - min(table[,hv],na.rm = TRUE))

        x <- lm( paste( hv , " ~ ",col_name," + " , paste( confounder_columns , collapse = " + ") ) , data = table , na.action = na.exclude)

        HV_SNP[i,1] <- summary(x)$coefficients[2,1]
        HV_SNP[i,2] <- max(table[,col_name] , na.rm=TRUE) - min(table[,col_name] , na.rm=TRUE)
        HV_SNP[i,3] <- summary(x)$coefficients[2,4]
        HV_SNP[i,4] <- summary(x)$r.squared
        i<-i+1
      }
      
      #print(paste(GENDER , HEMISPHERE))
      #print(HV_SNP)
      
#      if (GENDER=="Both Genders")
#        if (HEMISPHERE=="Both Hemispheres"){

      bplot <- barplot(main = paste(GENDER , HEMISPHERE) , xpd = FALSE ,ylim = c(summary(HV_SNP[,4])[1]-0.001,summary(HV_SNP[,4])[6]+0.1),
                       HV_SNP[,4] , xlab = "P-Value Threshold" , ylab = "R-Squared" , 
                       names=thresholds )
      
      text( x=c(1:14)*1.2,  y=HV_SNP[,4]+0.0001, formatC(HV_SNP[,3], format = "e", digits = 2) , srt=80 , adj=c(0,-1))
#        }
      max_col_name <-  names(HV_SNP[,4])[match(max(HV_SNP[,4]),HV_SNP[,4])]
      table <- table[  order(table[,max_col_name]) , ]
      hv <- table[ , hv]

#      if (GENDER=="Both Genders")
#        if (HEMISPHERE=="Both Hemispheres"){
          means <- 1:100
          percs <- 1:100
          for ( i in 1:100 ){
            # THIS METHOD CAlCULATESS MEANS AS PERCENTILE OF SAMPLES
            start_ind = nrow(table) * ((i-1)/100)
            end_ind = nrow(table) * (i/100)
            window <- hv[start_ind:end_ind]
            means[i] <- mean(window , na.rm = TRUE)

          }
          
          
#          var_vals <- 1:45
#          var_diff_old <- 0
#          for( i in 1:45 ){
#            means_without <- tail( head(means , length(means)-i) , length(means)-(2*i))
#            var_diff <- ( (var(means) - var(means_without)) / var(means) ) * 100
#            var_vals[i] <- var_diff-var_diff_old
#            print(paste("at",i,"incremntal var diff is",var_vals[i],"%"))
#            var_diff_old <- var_diff
#          }
          
          plot(means , main = paste("Percentiles for Threshold: " , max_col_name) )
          lines( predict(lm( head(means,-3) ~ poly(head(percs,-3),3) )) , col="grey" )
          
          
#          means <- means[!is.na(means)]
#          percent_of_range <- seq(0,50)*0
#          for( i in seq(1,45) ){
#          top_range <- range(head(means,i))[2] - range(head(means,i))[1]
#          bot_range <- range(tail(means,i))[2] - range(tail(means,i))[1]
#          full_range <- range(means)[2] - range(means)[1]
#          percent_of_range[i] <- (top_range + bot_range) / full_range
#          }
#          par(mfrow=c(1,1))
#          plot(percent_of_range , main = "PERCENT OF VARINACE EXPLAINED")
#        }
    }
  }
  return(HV_SNP)
}

