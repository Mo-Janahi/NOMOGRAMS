NOMOGRAMS

Project Goal / Idea :
	To build better Nomgrams by taking genetic information into account.

Steps in the script:

 1- Read in the UKBB table of 500k subjects
 2- Exclude those without imaging 
 3- Exclude those who have a history of head trauma or substance abuse 

 Now looking at Hippocampal volume:
 4- Exclude outliers (more than 5 MAE's away from mean). (Do this twice. Once for each hemisphere)
 5- Regress out the effect of: Head Scaling and Date of MRI Scan (Do this for each hemisphere seperately)

Using these corrected volumes:
 6- Generate the nomogram. (Sliding quantile window analysis to produce percentile curves)

 
For PRS generation/analysis:
 1- Create a .pheno file with 3 colmns: IID and FID are the "eid"s from the ukb table and 
    a pheno column which contains the sum of the left and right HV (corrected)

 2- Using that file, and the UKBB plink files, and the ENEGMA GWAS files, 
    input them into PRSsice to find some Polygenic Risk Scores (PRS)

 3- Back in R, read in the files and find the association between HV and the PRS scores found at
    Different thresholds with a linear regression model. Bar plot all results.

