rm(list = ls())
options(stringsAsFactors = FALSE)


library(data.table)
library(dplyr)


##for QC source step 1 file
source("step1_IQR_based_outlier_detection.R")
source("step2_GMM_based_biomarker_cutoff.R")
source("step3_SCeNS_calculation.R")


names(Input_file)

##for QC parameters needed are 
#dat="Your input file"
#x=biomarker column name, QC will be done column wise one biomarker per call to the function
##cols_keep= columns that you want to keep from you file, can also include your bimarker column; however the code keeps the biomarker column specified by "x" regardless
##and unique identifier to merge if mutiple markers are being QCed 
##log transformation is needed for this QC, this option is user defined and would be "TRUE" normally
##however if your data is in log scale already (eg Olink NPX or Alamar NPQ values) then use "FALSE"


Abeta42_check <- QC_function( dat= Input_file, ##input file to be QCed 
                       x="Abeta42", ##ypur biomarker value column name
                       cols_keep = c("SampleID", "Status"), ##which column to keep
                       apply_log = TRUE) #False = dont apply log ; here dat are synthetic Alamar NPQ values so already in log2 scale



pTau217_check<- QC_function( dat= Input_file, ##input file to be QCed 
                             x="pTau217", ##ypur biomarker value column name
                             cols_keep = c("SampleID", "Status"), ##which column to keep
                             apply_log = TRUE) #False = dont apply log ; here dat are synthetic Alamar NPQ values so already in log2 scale

pTau181_check<- QC_function( dat= Input_file, ##input file to be QCed 
                             x="pTau181", ##ypur biomarker value column name
                             cols_keep = c("SampleID", "Status"), ##which column to keep
                             apply_log = TRUE) #False = dont apply log ; here dat are synthetic Alamar NPQ values so already in log2 scale



###################now perform cutoff generation############################
##please note that because of disimilar direction of change of ABETA42 and ptau species in fluid samples, seperate function for each

##users need to provide data file name, raw column name and z score name


##this is done post QC 
names(Abeta42_check) 


Get_ABeta_cutoff<- perform_AB_Dichotomization(x=Abeta42_check, 
                                              raw_col="Abeta42", ##QCed raw value column 
                                              zscore_col="Abeta42_Zscore") ##z score of abeta42 calculated following QC in above steps

##this object contains the dataframe with biomarker status label
##plot showing the identified distribution 
##z score biomarker cutoff ## also printed out in console

##access each component as

Get_ABeta_cutoff_dataframe<- Get_ABeta_cutoff[[1]]
Get_ABeta_cutoff_plot<- Get_ABeta_cutoff[[2]]
Abeta_cutoff<- Get_ABeta_cutoff[[3]] ## this will be needed for our SCeNS calculation


##for ptau217
names(pTau217_check)

Get_ptau217_cutoff<- perform_pTau_Dichotomization(x=pTau217_check, 
                                              raw_col= "pTau217", ##QCed raw value column 
                                              zscore_col="pTau217_Zscore") ##z score of abeta42 calculated following QC in above steps

##this object contains the dataframe with biomarker status label
##plot showing the identified distribution 

##access each component as

Get_ptau217_cutoff_dataframe<- Get_ptau217_cutoff[[1]]
Get_ptau217_cutoff_plot<-Get_ptau217_cutoff[[2]]
pTau217_cutoff<- Get_ptau217_cutoff[[3]]

##pTau181

names(pTau181_check)

Get_ptau181_cutoff<- perform_pTau_Dichotomization(x=pTau181_check, 
                                                  raw_col= "pTau181", ##QCed raw value column 
                                                  zscore_col="pTau181_Zscore") ##z score of abeta42 calculated following QC in above steps

##this object contains the dataframe with biomarker status label
##plot showing the identified distribution 

##access each component as

Get_ptau181_cutoff_dataframe<- Get_ptau181_cutoff[[1]]
Get_ptau181_cutoff_plot<- Get_ptau181_cutoff[[2]]
pTau181_cutoff<- Get_ptau181_cutoff[[3]]



Get_ABeta_cutoff_plot<-  Get_ABeta_cutoff_plot +
  annotate("text", x = 0, y = 0.6, 
           label = "Gaussian Mixture Density identified in AB42",
           fontface="bold")


Get_ptau181_cutoff_plot<-  Get_ptau181_cutoff_plot +
  annotate("text", x = 0, y = 0.5, 
           label = "Gaussian Mixture Density identified in pTau-181",
           fontface="bold")


Get_ptau217_cutoff_plot <-  Get_ptau217_cutoff_plot +
  annotate("text", x = 0, y = 0.5, 
           label = "Gaussian Mixture Density identified in pTau-217",
           fontface="bold")


#############Now SCeNS###################
##three parameter that need to be supplied are dataframe to use, column to use and cutoff to center the values (As determined by GMM above)
##SCeNS is designe dto work with zscores and z-score based cutoff supply the z score column as column to use


AB42_SCeNS_calculation<- AB42_recalib(data=Get_ABeta_cutoff_dataframe, 
                                      column_to_use= "Abeta42_Zscore", 
                                      cutoff=Abeta_cutoff) ##stored above from GMM run


pTau_SCeNS_calculation<- pTau_recalib(data=Get_ptau181_cutoff_dataframe, 
                                      column_to_use= "pTau181_Zscore", 
                                      cutoff=pTau181_cutoff) ##stored above from GMM run


pltAB_1<- ggplot(data= AB42_SCeNS_calculation, aes(x= as.numeric(SCeNS_AB42), fill = Abeta42_Group))+
  geom_histogram(col="white")+
  labs(title="Distribution of re-scaled SCeNS", x= "AB42 SCeNS score", y= "Biomarker status") +
  scale_fill_manual(values = c( "steelblue","salmon"))+
  theme_classic()+
  guides(fill=guide_legend(title="Status"))


pltAB_2<- ggplot(data= AB42_SCeNS_calculation, aes(x= Abeta42_Group,y= as.numeric(SCeNS_AB42),fill = Abeta42_Group))+
  geom_boxplot()+
  labs(title ="SCeNS score seperation between biomarker status group",x= "Biomarker status", y= "AB42 SCeNS score")+
  geom_hline(yintercept = 0, lty=2, color="black")+
  scale_fill_manual(values = c( "steelblue","salmon"))+
  theme_classic()+
  guides(fill=guide_legend(title="Status"))

pltAB_3<- ggplot(data= AB42_SCeNS_calculation, aes(x= as.numeric(SCeNS_AB42), y =as.numeric(Abeta42_Zscore), fill = Abeta42_Group, color=Abeta42_Group))+
  geom_point(shape=21)+
  labs(title="Correlation between SCeNS and underlying Z-score", x= "SCeNS score", y= "Z scores")+
  scale_fill_manual(values = c( "steelblue","salmon"))+
  scale_color_manual(values = c( "steelblue","salmon"), guide = "none")+
  theme_classic()+
  guides(fill=guide_legend(title="Status"))

