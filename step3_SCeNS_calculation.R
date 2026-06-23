
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)


AB42_recalib<- function(data, column_to_use, cutoff){
  data <- as.data.frame(data)
  data<-  data %>% filter(!is.na(.data[[column_to_use]]))
  ##recalibrtae the range of distribution
  data <- data %>% mutate(Recalib_AB42 = .data[[column_to_use]] - cutoff)
  scale_factor <- (-1) * 100 / max(abs(data$Recalib_AB42))  ##-1 because lower z score mean AB+ and higher zscore mean AB-
  ##rescale
  data$SCeNS_AB42 <- scale_factor * data$Recalib_AB42
  data$Recalib_AB42<- NULL 
  return(data)
  }


pTau_recalib<- function(data, column_to_use, cutoff){
  data= data %>% filter(!is.na(.data[[column_to_use]]))
  data <- data %>% 
    mutate(Recalib_pTau = .data[[column_to_use]] - cutoff)
  ##recalibrtae the range of distribution 
  scale_factor<- 1*100/max(abs(data$Recalib_pTau)) 
  ##rescale
  data$Standard_pTau <- scale_factor * data$Recalib_pTau
  data$Recalib_pTau <- NULL
  return(data)
}

