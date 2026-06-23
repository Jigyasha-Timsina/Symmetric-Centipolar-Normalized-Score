options(stringsAsFactors = FALSE)



library(data.table)
library(dplyr)
library(mclust)


perform_AB_Dichotomization <- function(x, raw_col, zscore_col) {
  
  keep_cols<- names(x) ##later to keep required cols
  
  ### GMM clustering into two groups using z-score
  set.seed(1)
  lc_out <- Mclust(x[[zscore_col]], G = 2)
  print(summary(lc_out))
  
  ## Add grouping into main data table
  x <- data.frame(x, lc_out$z, lc_out$classification)
  
  ##find which value is cutoff
  per50_df <- x[which.min(abs(0.50 - x$X1)), ]
  cutoff_zscore <- per50_df[[zscore_col]]
  cutoff_raw <- per50_df[[raw_col]]
  print(paste0("AB42 z-score cut-off :", cutoff_zscore))
  print(paste0("AB42 rawvalue cut-off :", cutoff_raw))
  
  # Classify based on z-score cutoff
  x$Amyloid_bin <- NA
  x$Amyloid_bin[x[[zscore_col]] <  cutoff_zscore] <- 1  # Amyloid+
  x$Amyloid_bin[x[[zscore_col]] >= cutoff_zscore] <- 2  # Amyloid-
  print(paste0("A+ subjects : ", table(x$Amyloid_bin)[1]))
  print(paste0("A- subjects : ", table(x$Amyloid_bin)[2]))
  
  ##assign A=/- label
  group_col <- paste0(raw_col, "_Group")
  x[[group_col]] <- ifelse(x$Amyloid_bin == 1, "A+", "A-")
  
  x<- x %>% dplyr::select(all_of(c(keep_cols,group_col)))
  
  ### Density overlay parameters
  val  <- length(lc_out$param$variance$sigmasq)
  iter <- seq(-3, 3, 0.5)
  y1   <- lc_out$param$pro[1] * dnorm(iter, lc_out$param$mean[1], sqrt(lc_out$param$variance$sigmasq[1]))
  y2   <- lc_out$param$pro[2] * dnorm(iter, lc_out$param$mean[2], sqrt(lc_out$param$variance$sigmasq[val]))
  
  density_line1 <- data.frame(iter, y1)
  density_line2 <- data.frame(iter, y2)
  
  # Temp column for ggplot aes()
  x$.zscore_plot <- x[[zscore_col]]
  
  plot_distribution <- ggplot() +
    geom_histogram(data  = x, aes(x = .zscore_plot, y = ..density..), color = "#440154ff", fill = "white") +
    geom_vline(data = x, aes(xintercept = mean(.zscore_plot)),lty = "dashed") +
    geom_line(data = density_line1, aes(x = iter, y = y1),
              color = "red",  lty = "dashed", lwd = 0.9) +
    geom_line(data = density_line2, aes(x = iter, y = y2),
              color = "blue", lty = "dashed", lwd = 0.9) +
    labs(x = paste0(raw_col, " (z-score)")) +
    theme_classic()
  
  # Remove temp column before returning
  x$.zscore_plot <- NULL
  
  output <- list(x, plot_distribution,cutoff_zscore)
  return(output)
  }


####################################becasue direction differs between the two###################
##seperate fucntion for tau and ptau species

perform_pTau_Dichotomization <- function(x, raw_col, zscore_col) {
  keep_cols<- names(x) ##later to keep required cols
  
  ### GMM clustering into two groups using z-score
  set.seed(1)
  lc_out <- Mclust(x[[zscore_col]], G = 2)
  print(summary(lc_out))
  
  ## Add grouping into main data table
  x <- data.frame(x, lc_out$z, lc_out$classification)
  
  ##find which value is cutoff
  per50_df      <- x[which.min(abs(0.50 - x$X1)), ]
  cutoff_zscore <- per50_df[[zscore_col]]
  cutoff_raw    <- per50_df[[raw_col]]
  print(paste0("pTau z-score cut-off :", cutoff_zscore))
  print(paste0("pTau rawvalue cut-off :", cutoff_raw))
  
  # Classify based on z-score cutoff
  x$pTau_bin <- NA
  x$pTau_bin[x[[zscore_col]] <  cutoff_zscore] <- 1  # pTau-
  x$pTau_bin[x[[zscore_col]] >= cutoff_zscore] <- 2  # pTau+
  print(paste0("T- subjects : ", table(x$pTau_bin)[1]))
  print(paste0("T+ subjects : ", table(x$pTau_bin)[2]))
  
  ##assign T+/- label
  group_col <- paste0(raw_col, "_Group")
  x[[group_col]] <- ifelse(x$pTau_bin == 1, "T-", "T+")
  
  x<- x %>% dplyr::select(all_of(c(keep_cols,group_col)))
  
  ### Density overlay parameters
  val  <- length(lc_out$param$variance$sigmasq)
  iter <- seq(-3, 3, 0.5)
  y1   <- lc_out$param$pro[1] * dnorm(iter, lc_out$param$mean[1], sqrt(lc_out$param$variance$sigmasq[1]))
  y2   <- lc_out$param$pro[2] * dnorm(iter, lc_out$param$mean[2], sqrt(lc_out$param$variance$sigmasq[val]))
  
  density_line1 <- data.frame(iter, y1)
  density_line2 <- data.frame(iter, y2)
  
  # Temp column for ggplot aes()
  x$.zscore_plot <- x[[zscore_col]]
  
  plot_distribution <- ggplot() +
    geom_histogram(data  = x, aes(x = .zscore_plot, y = ..density..), color = "#440154ff", fill = "white") +
    geom_vline(data = x, aes(xintercept = mean(.zscore_plot)), lty = "dashed") +
    geom_line(data = density_line1, aes(x = iter, y = y1),
              color = "blue",  lty = "dashed", lwd = 0.9) +
    geom_line(data = density_line2, aes(x = iter, y = y2),
              color = "red", lty = "dashed", lwd = 0.9) +
    labs(x = paste0(raw_col, " (z-score)")) +
    ylim(0, 0.5) +
    theme_classic()
  
  # Remove temp column before returning
  x$.zscore_plot <- NULL
  
  output <- list(x, plot_distribution,cutoff_zscore)
  return(output)
}
    
