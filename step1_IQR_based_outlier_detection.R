options(stringsAsFactors = FALSE)


library(data.table)
library(dplyr)
library(lubridate)
library(readxl)
library(stringi)
library(stringr)
library(tidyr)


##QC Function: this fucntion performs optional log10, IQR outlier removal followed by  Z-score normalization
##returns: data frame with outliers removed and a new Z-score column appended.


##info needed from users
#dat="Your input file name"
#x=biomarker column name, QC will be done column wise one biomarker per call to the function
##cols_keep= columns that you want to keep from you file, can also include your bimarker column; however the code keeps the biomarker column specified by "x" regardless
##and unique identifier to merge if mutiple markers are being QCed 
##log transformation is needed for this QC, this option is user defined and would be "TRUE" normally
##however if your data is in log scale already (eg Olink NPX or Alamar NPQ values) then use "FALSE"




##wrapper for QC function

QC_function<- function(dat,x, cols_keep, apply_log, cut_off = 1.5) {
  
  dat <- as.data.frame(dat)
  ##make sure biomarker column is present in the matrix and is numeric values
  x <- as.character(x)
  if (!x %in% names(dat)) stop(paste0("Column '", x, "' not found in data."))
  if (!is.numeric(dat[[x]])) stop(paste0("Column '", x, "' must be numeric."))
  
  ##chekc if any column asked to be kept not in the dataframe and warn the user
  if (!is.null(cols_keep)) {
    missing_cols <- setdiff(cols_keep, names(dat))
    if (length(missing_cols) > 0)
      warning(paste0("These cols_keep columns were not found and will be ignored: ", paste(missing_cols, collapse = ", ")))
    cols_keep <- union(intersect(cols_keep, names(dat)), x)
    dat <- dat[, cols_keep, drop = FALSE]}
  
  ##check if there are NAs in the matrix and remove if they do
  print(paste0("Missing ",x," :", sum(is.na(dat[[x]]))))
  dat <- tidyr::drop_na(dat, all_of(x))
  
  ##this step is for log transformation ; if apply_log=FALSE log wont be applied
  ##by default we use log 10 if needed change the log base in code
  if (apply_log) {
    if (any(dat[[x]] <= 0, na.rm = TRUE))
      warning("Non-positive values detected; log10 of these will be NaN/Inf.")
    dat[[paste0(x, "_log10")]] <- log10(dat[[x]])
    work_col_name <- paste0(x, "_log10")
  } else {
    work_col_name <- x
  }
  
  ##check for outlier values based on IQR using the log scale values and remove outlier values
  dat$outlier_status <- NA
  IQR <- IQR(dat[[work_col_name]])
  Quantiles <- quantile(dat[[work_col_name]])
  
  for (i in 1:nrow(dat)) {
    if (dat[[work_col_name]][i] > (Quantiles['75%'] + IQR * cut_off)) { 
      dat$outlier_status[i] <- "Outlier"
    } else if (dat[[work_col_name]][i] < (Quantiles['25%'] - IQR * cut_off)) { 
        dat$outlier_status[i] <- "Outlier"
    } else { 
          dat$outlier_status[i] <- "Not Outlier"
          }
    }
  
  print(paste0("No. of Outliers: ", table(dat$outlier_status)["Outlier"]))
  dat <- dat %>% filter(!outlier_status == "Outlier")
  dat$outlier_status <- NULL
  
  ##Normalize the values
  zscore_col <- paste0(x,"_Zscore")
  dat[[zscore_col]] <- as.numeric(scale(dat[[work_col_name]], center = TRUE, scale = TRUE))
  return(dat)
  }




