library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dtw)
library(dbscan)
library(FNN)
library(caTools)
library(Metrics)
library(rsample)
library(Rfast)
library(lubridate)
library(patchwork)
library(tibble)
library(kernlab)
library(dtwclust)
library(Rcpp)
library(Metrics)
library(dplyr)
library(tidyverse)
library(scoringRules)
library(foreach)
library(doParallel)
library(kernlab)
library(Matrix)
library(MASS)
library(progress)

####################### KERNELS 

data <- read.table("DATA.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meteo  <- read.table("METEO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geno   <- read.table("GENO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

data$yield <- (data$yield - mean(data$yield))/sd(data$yield)
data$prot <- (data$prot - mean(data$prot))/sd(data$prot)

# normalize the different weather variables
normalize_zscore <- function(x) {
  y <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) 
  if(sd(x, na.rm=TRUE)==0){y <- 0}
  return(y)}
normalize_min_max <- function(x) {
  y <- (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  if(max(x, na.rm=TRUE)==0 & min(x, na.rm=TRUE)){y <- 0}
  return(y)}

c2 <- ncol(meteo)
meteo[,-c2] <- lapply(meteo[,-c2], normalize_zscore)
varNames <- sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()
meteo1 <- meteo %>%
  dplyr::select(matches(paste(varNames[1:7], collapse = "|")))
meteo <- cbind(meteo$Env, meteo1)
colnames(meteo)[1] <- 'Env'

#################################################################################
# HAMMING for GENO

# count the number of equal or NA entries
hamming_distance_matrix <- function(maat) {
  n <- nrow(maat)
  hamming_matrix <- outer(1:n, 1:n, Vectorize(function(i, j) {
    v <- (maat[i, ] == maat[j, ]) | is.na(maat[i, ]) | is.na(maat[j, ]) 
    return((ncol(maat) - sum(v)) / ncol(maat)) 
  }))
  return(hamming_matrix)
}
# geno_dist_hamming <- hamming_distance_matrix(geno[,-1])
# rownames(geno_dist_hamming) <- geno$variety_name
# colnames(geno_dist_hamming) <- geno$variety_name
# write.table(geno_dist_hamming, file = "GENO_HAMMING.txt", sep = "\t", row.names = FALSE, quote = FALSE)
geno_dist_hamming   <- read.table("GENO_HAMMING.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)


#################################################################################
# Alignment for meteo  (from Marco Cuturi's webpage)
# Useful constants
LOG0 <- -10000  # log(0)

# LOGP function: stable computation of log(exp(x) + exp(y))
LOGP <- function(x, y) {
  if (x > y) {
    return(x + log1p(exp(y - x)))
  } else {
    return(y + log1p(exp(x - y)))
  }
}

# Global Alignment Kernel function
logGAK <- function(seq1, seq2, sigma, triangular = 0) {
  nX <- nrow(seq1)  # Number of rows in seq1
  nY <- nrow(seq2)  # Number of rows in seq2
  dimvect <- ncol(seq1)  # Dimension of the time series vectors
  
  # Length of a column for dynamic programming
  cl <- nY + 1
  
  # Initialize logM to store two successive columns (dynamic programming table)
  logM <- matrix(LOG0, nrow = cl, ncol = 2)
  logM[1, 1] <- 0  # Initialize the lower-left cell with log(1) = 0
  
  # Maximum of abs(i - j) when 1 <= i <= nX and 1 <= j <= nY
  trimax <- max(nX - 1, nY - 1)
  
  # Triangular coefficients initialization
  logTriangularCoefficients <- rep(0, trimax + 1)
  if (triangular > 0) {
    for (i in seq_len(min(trimax + 1, triangular)) - 1) {
      logTriangularCoefficients[i + 1] <- log(1 - i / triangular)
    }
  }
  
  # Sigma factor for Gaussian kernel
  Sig <- -1 / (2 * sigma^2)
  
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i, ] - seq2[j, ])^2)
        gram <- logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
        gram <- gram - log(2 - exp(gram))
        
        # Update logM for dynamic programming
        frompos1 <- logM[j + 1, old]
        frompos2 <- logM[j, cur]
        frompos3 <- logM[j, old]
        aux <- LOGP(frompos1, frompos2)
        logM[j + 1, cur] <- LOGP(aux, frompos3) + gram
      }
    }
    # Swap cur and old
    cur <- 3 - cur
    old <- 3 - old
  }
  
  # Return the final result
  return(logM[nY + 1, old])
}

logGAK_derivative <- function(seq1, seq2, sigma, triangular = 0) {
  nX <- nrow(seq1)  # Number of rows in seq1
  nY <- nrow(seq2)  # Number of rows in seq2
  dimvect <- ncol(seq1)  # Dimension of the time series vectors
  
  # Length of a column for dynamic programming
  cl <- nY + 1
  
  # Initialize logM to store two successive columns (dynamic programming table)
  logM <- matrix(LOG0, nrow = cl, ncol = 2)
  logM[1, 1] <- 0  # Initialize the lower-left cell with log(1) = 0
  
  # Maximum of abs(i - j) when 1 <= i <= nX and 1 <= j <= nY
  trimax <- max(nX - 1, nY - 1)
  
  # Triangular coefficients initialization
  logTriangularCoefficients <- rep(0, trimax + 1)
  if (triangular > 0) {
    for (i in seq_len(min(trimax + 1, triangular)) - 1) {
      logTriangularCoefficients[i + 1] <- log(1 - i / triangular)
    }
  }
  
  # Sigma factor for Gaussian kernel
  Sig <- -1 / (2 * sigma^2)
  del_Sig=sigma^(-3)
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i, ] - seq2[j, ])^2)
        gram <- logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
        gram <- gram - log(2 - exp(gram))
        del_gram=diff_sq * del_Sig*(1+1/(2-exp(gram))*exp(gram))
        
        # Update logM for dynamic programming
        frompos1 <- logM[j + 1, old]
        frompos2 <- logM[j, cur]
        frompos3 <- logM[j, old]
        aux <- LOGP(frompos1, frompos2)
        logM[j + 1, cur] <- LOGP(aux, frompos3) + del_gram
      }
    }
    # Swap cur and old
    cur <- 3 - cur
    old <- 3 - old
  }
  
  # Return the final result
  return(logM[nY + 1, old])
}

varNames <- sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()


################################################################################


# setupp
keri <- 1  # which kernel combination
setupp <- '1' # new variety or new environment
leakage <- 'no'
leak_spec <- 'random'
noisy_obs = TRUE
train_prop = 0.8
nr_outer_split = 1
lr=.01
max_iter=20
notcluster=1
# setwd("C:/Users/leafr/OneDrive - Universitaet Bern/Desktop/Projects/Agroscope/1_DEF")
data <- read.table("DATA.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meteo  <- read.table("METEO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geno   <- read.table("GENO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# normalize the different weather variables
normalize_zscore <- function(x) {
  y <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) 
  if(sd(x, na.rm=TRUE)==0){y <- 0}
  return(y)}
normalize_min_max <- function(x) {
  y <- (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  if(max(x, na.rm=TRUE)==0 & min(x, na.rm=TRUE)){y <- 0}
  return(y)}

c2 <- ncol(meteo)
meteo[,-c2] <- lapply(meteo[,-c2], normalize_zscore)
varNames <- sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()
meteo1 <- meteo[,which(grepl(paste(varNames[1:7], collapse = "|"), colnames(meteo)))]
meteo <- cbind(meteo$Env, meteo1)
colnames(meteo)[1] <- 'Env'

###################################################################################
# create the kernel matrices

create_kernelmat <- function(kernel_meteo, kernel_geno, i1, i2, mat){
  if(kernel_meteo == 'exp'){
    meteo_dist_eucl <- as.matrix(dist(meteo[,-1], method = "euclidean"))
    rownames(meteo_dist_eucl) <- meteo$Env
    colnames(meteo_dist_eucl) <- meteo$Env
    int <- seq(0.1, 1.9, length.out = 10)[i1]*max(meteo_dist_eucl)   # c(0.2,0.3,0.4,0.5,0.6)[i1]
    Dmat <- meteo_dist_eucl[mat$Env, mat$Env]
    Kmat_meteo <- exp(-Dmat/(int))}
  
  # alignments kernel
  if(kernel_meteo == 'ali'){
    s <- i1 # seq(2,8,length.out=10)[i1]
    meteo_ali1   <- read.table(paste0("METEO_GA_sigma__", s, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    rownames(meteo_ali1) <- colnames(meteo_ali1)
    meteo_ali <- meteo_ali1[mat$Env, mat$Env]
    rownames(meteo_ali) <- colnames(meteo_ali)  
    Kmat_meteo <- meteo_ali}
  
  if(kernel_geno == 'ham'){
    geno_dist_hamming   <- read.table("GENO_HAMMING.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)
    int <- seq(0.1, 1.9, length.out = 10)[i2]*max(geno_dist_hamming)   # seq(0.1, 1.5, length.out = 10)[i2]*max(geno_dist_hamming)
    Dmat <- geno_dist_hamming[mat$variety_name, mat$variety_name]
    Kmat_geno <- exp(-Dmat/(int))}
  
  # spectrum kernel
  if(kernel_geno == 'spe'){
    s <- i2
    geno_spec1   <- read.table(paste0("GENO_SPEC", s, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_spec1)[colnames(geno_spec1) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_spec1)[colnames(geno_spec1) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_spec1) <- colnames(geno_spec1)
    geno_spec <- geno_spec1[mat$variety_name, mat$variety_name]/max(geno_spec1)
    rownames(geno_spec) <- colnames(geno_spec)
    Kmat_geno <- geno_spec}
  
  Kmat <- Kmat_meteo*Kmat_geno
  # Kmat <- Kmat + 1e-3*diag(1,nrow(mat))
  return(Kmat)}

create_kernelmat_continuous <- function(kernel_meteo, kernel_geno, theta1,theta2,k,
                                        mat){
  if(kernel_meteo == 'exp'){
    meteo_dist_eucl <- as.matrix(dist(meteo[,-1], method = "euclidean"))
    rownames(meteo_dist_eucl) <- meteo$Env
    colnames(meteo_dist_eucl) <- meteo$Env
    #int <- seq(0.1, 1.9, length.out = 10)[i1]*max(meteo_dist_eucl)   # c(0.2,0.3,0.4,0.5,0.6)[i1]
    int <- theta1*max(meteo_dist_eucl)
    Dmat <- meteo_dist_eucl[mat$Env, mat$Env]
    Kmat_meteo <- exp(-Dmat/(int))
  }
  
  # alignments kernel
  if(kernel_meteo == 'ali'){
    # s <- i1 # seq(2,8,length.out=10)[i1
    #dat=meteo[ind,]
    Kernel_GA <- matrix(NA, nrow=nrow(meteo), ncol=nrow(meteo))
    
    for(j1 in 1:nrow(meteo)){
      for(j2 in 1:j1){
        
        meteoM1 <- matrix(NA, nrow = 6, ncol = 7)
        meteoM2 <- matrix(NA, nrow = 6, ncol = 7)
        for(jj in 1:7){
          mm1 <- as.numeric(meteo[j1,grep(varNames[jj], colnames(meteo))])
          mm2 <- as.numeric(meteo[j2,grep(varNames[jj], colnames(meteo))])
          meteoM1[1:length(mm1),jj] <- mm1
          meteoM2[1:length(mm2),jj] <- mm2
        }
        
        t1 <- logGAK(meteoM1, meteoM2, sigma = theta1, triangular = 0)
        t2 <- logGAK(meteoM2, meteoM2, sigma = theta1, triangular = 0)
        t3 <- logGAK(meteoM1, meteoM1, sigma = theta1, triangular= 0)
        
        Kernel_GA[j1,j2] <- exp(t1-.5*(t2+t3))
        Kernel_GA[j2,j1] <- exp(t1-.5*(t2+t3))
        
      }}
    
    rownames(Kernel_GA) <- meteo$Env
    colnames(Kernel_GA) <- meteo$Env
    
    Kmat_meteo <-Kernel_GA[mat$Env, mat$Env]}
  
  if(kernel_geno == 'ham'){
    geno_dist_hamming   <- read.table("GENO_HAMMING.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)
    int <- theta2*max(geno_dist_hamming)   # seq(0.1, 1.5, length.out = 10)[i2]*max(geno_dist_hamming)
    Dmat <- geno_dist_hamming[mat$variety_name, mat$variety_name]
    Kmat_geno <- exp(-Dmat/(int))}
  
  # spectrum kernel
  if(kernel_geno == 'spe'){
    s <- k
    geno_spec1   <- read.table(paste0("GENO_SPEC", s, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_spec1)[colnames(geno_spec1) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_spec1)[colnames(geno_spec1) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_spec1) <- colnames(geno_spec1)
    geno_spec <- geno_spec1[mat$variety_name, mat$variety_name]/max(geno_spec1)
    rownames(geno_spec) <- colnames(geno_spec)
    Kmat_geno <- geno_spec}
  
  Kmat <- Kmat_meteo*Kmat_geno
  # Kmat <- Kmat + 1e-3*diag(1,nrow(mat))
  return(Kmat)}


create_kernelmat_derivative<- function(kernel_meteo, kernel_geno, theta1,theta2, k,mat){
  if(kernel_meteo == 'exp'){
    meteo_dist_eucl <- as.matrix(dist(meteo[,-1], method = "euclidean"))
    rownames(meteo_dist_eucl) <- meteo$Env
    colnames(meteo_dist_eucl) <- meteo$Env
    int <- theta1*max(meteo_dist_eucl)   # c(0.2,0.3,0.4,0.5,0.6)[i1]
    Dmat <- meteo_dist_eucl[mat$Env, mat$Env]
    Kmat_meteo <- exp(-Dmat/(int))
    Kmat_meteo_deriv <- exp(-Dmat/(int))*Dmat/(int^2)*max(meteo_dist_eucl)}
  
  # alignments kernel
  if(kernel_meteo == 'ali'){
    
    Kernel_GA=Kernel_GA_deriv <- matrix(NA, nrow=nrow(meteo), ncol=nrow(meteo))
    rownames(Kernel_GA) <- meteo$Env
    colnames(Kernel_GA) <- meteo$Env
    rownames(Kernel_GA_deriv) <- meteo$Env
    colnames(Kernel_GA_deriv) <- meteo$Env
    
    for(j1 in 1:nrow(meteo)){
      for(j2 in 1:j1){
        
        meteoM1 <- matrix(NA, nrow = 6, ncol = 7)
        meteoM2 <- matrix(NA, nrow = 6, ncol = 7)
        for(jj in 1:7){
          mm1 <- as.numeric(meteo[j1,grep(varNames[jj], colnames(meteo))])
          mm2 <- as.numeric(meteo[j2,grep(varNames[jj], colnames(meteo))])
          meteoM1[1:length(mm1),jj] <- mm1
          meteoM2[1:length(mm2),jj] <- mm2
        }
        
        t1 <- logGAK(meteoM1, meteoM2, sigma = theta1, triangular = 0)
        t1_deriv=logGAK_derivative(meteoM1, meteoM2, sigma = theta1, triangular = 0)
        t2 <- logGAK(meteoM2, meteoM2, sigma = theta1, triangular = 0)
        t2_deriv=logGAK_derivative(meteoM2, meteoM2, sigma = theta1, triangular = 0)
        
        t3 <- logGAK(meteoM1, meteoM1, sigma = theta1, triangular= 0)
        t3_deriv=logGAK_derivative(meteoM1, meteoM1, sigma = theta1, triangular = 0)
        
        
        Kernel_GA_deriv[j1,j2] <- exp(t1-.5*(t2+t3))*(t1_deriv-.5*(t2_deriv+t3_deriv))
        
        Kernel_GA_deriv[j2,j1] <- exp(t1-.5*(t2+t3))*(t1_deriv-.5*(t2_deriv+t3_deriv))
        Kernel_GA[j1,j2] <- exp(t1-.5*(t2+t3))*(t1_deriv-.5*(t2_deriv+t3_deriv))
        
        Kernel_GA[j2,j1] <- exp(t1-.5*(t2+t3))
        
      }}
    
    
    Kmat_meteo_deriv <-Kernel_GA_deriv[mat$Env, mat$Env]
    Kmat_meteo <-Kernel_GA[mat$Env, mat$Env]
  }
  
  if(kernel_geno == 'ham'){
    geno_dist_hamming   <- read.table("GENO_HAMMING.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)
    int <- theta2*max(geno_dist_hamming)   # seq(0.1, 1.5, length.out = 10)[i2]*max(geno_dist_hamming)
    Dmat <- geno_dist_hamming[mat$variety_name, mat$variety_name]
    Kmat_geno <- exp(-Dmat/(int))
    Kmat_geno_deriv <- Kmat_geno*Dmat/(int^2)*max(geno_dist_hamming) }
  
  # spectrum kernel
  if(kernel_geno == 'spe'){
    s <- k
    geno_spec1   <- read.table(paste0("GENO_SPEC", s, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_spec1)[colnames(geno_spec1) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_spec1)[colnames(geno_spec1) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_spec1) <- colnames(geno_spec1)
    geno_spec <- geno_spec1[mat$variety_name, mat$variety_name]/max(geno_spec1)
    rownames(geno_spec) <- colnames(geno_spec)
    Kmat_geno_deriv= Kmat_geno <- geno_spec
  }
  
  
  # Kmat <- Kmat + 1e-3*diag(1,nrow(mat))
  return(list(del_theta1=Kmat_meteo_deriv*Kmat_geno,
              del_theta2=Kmat_meteo*Kmat_geno_deriv))
}

create_kernelmat_derivative<- function(kernel_meteo, kernel_geno, theta1,theta2, k,mat){
  if(kernel_meteo == 'exp'){
    meteo_dist_eucl <- as.matrix(dist(meteo[,-1], method = "euclidean"))
    rownames(meteo_dist_eucl) <- meteo$Env
    colnames(meteo_dist_eucl) <- meteo$Env
    int <- theta1*max(meteo_dist_eucl)   # c(0.2,0.3,0.4,0.5,0.6)[i1]
    Dmat <- meteo_dist_eucl[mat$Env, mat$Env]
    Kmat_meteo <- exp(-Dmat/(int))
    Kmat_meteo_deriv <- exp(-Dmat/(int))*Dmat/(int^2)*max(meteo_dist_eucl)}
  
  # alignments kernel
  if(kernel_meteo == 'ali'){
    
    Kernel_GA=Kernel_GA_deriv <- matrix(NA, nrow=nrow(meteo), ncol=nrow(meteo))
    rownames(Kernel_GA) <- meteo$Env
    colnames(Kernel_GA) <- meteo$Env
    rownames(Kernel_GA_deriv) <- meteo$Env
    colnames(Kernel_GA_deriv) <- meteo$Env
    
    for(j1 in 1:nrow(meteo)){
      for(j2 in 1:j1){
        
        meteoM1 <- matrix(NA, nrow = 6, ncol = 7)
        meteoM2 <- matrix(NA, nrow = 6, ncol = 7)
        for(jj in 1:7){
          mm1 <- as.numeric(meteo[j1,grep(varNames[jj], colnames(meteo))])
          mm2 <- as.numeric(meteo[j2,grep(varNames[jj], colnames(meteo))])
          meteoM1[1:length(mm1),jj] <- mm1
          meteoM2[1:length(mm2),jj] <- mm2
        }
        
        t1 <- logGAK(meteoM1, meteoM2, sigma = theta1, triangular = 0)
        t1_deriv=logGAK_derivative(meteoM1, meteoM2, sigma = theta1, triangular = 0)
        t2 <- logGAK(meteoM2, meteoM2, sigma = theta1, triangular = 0)
        t2_deriv=logGAK_derivative(meteoM2, meteoM2, sigma = theta1, triangular = 0)
        
        t3 <- logGAK(meteoM1, meteoM1, sigma = theta1, triangular= 0)
        t3_deriv=logGAK_derivative(meteoM1, meteoM1, sigma = theta1, triangular = 0)
        
        
        Kernel_GA_deriv[j1,j2] <- exp(t1-.5*(t2+t3))*(t1_deriv-.5*(t2_deriv+t3_deriv))
        
        Kernel_GA_deriv[j2,j1] <- exp(t1-.5*(t2+t3))*(t1_deriv-.5*(t2_deriv+t3_deriv))
        Kernel_GA[j1,j2] <- exp(t1-.5*(t2+t3))*(t1_deriv-.5*(t2_deriv+t3_deriv))
        
        Kernel_GA[j2,j1] <- exp(t1-.5*(t2+t3))
        
      }}
    
    
    Kmat_meteo_deriv <-Kernel_GA_deriv[mat$Env, mat$Env]
    Kmat_meteo <-Kernel_GA[mat$Env, mat$Env]
  }
  
  if(kernel_geno == 'ham'){
    geno_dist_hamming   <- read.table("GENO_HAMMING.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)
    int <- theta2*max(geno_dist_hamming)   # seq(0.1, 1.5, length.out = 10)[i2]*max(geno_dist_hamming)
    Dmat <- geno_dist_hamming[mat$variety_name, mat$variety_name]
    Kmat_geno <- exp(-Dmat/(int))
    Kmat_geno_deriv <- Kmat_geno*Dmat/(int^2)*max(geno_dist_hamming) }
  
  # spectrum kernel
  if(kernel_geno == 'spe'){
    s <- k
    geno_spec1   <- read.table(paste0("GENO_SPEC", s, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    colnames(geno_spec1)[colnames(geno_spec1) == "SCARO_.ex_SARO."] <- "SCARO_(ex_SARO)"
    colnames(geno_spec1)[colnames(geno_spec1) == "EMBLEM_.NEW."] <- "EMBLEM_(NEW)"
    rownames(geno_spec1) <- colnames(geno_spec1)
    geno_spec <- geno_spec1[mat$variety_name, mat$variety_name]/max(geno_spec1)
    rownames(geno_spec) <- colnames(geno_spec)
    Kmat_geno_deriv= Kmat_geno <- geno_spec
  }
  
  
  # Kmat <- Kmat + 1e-3*diag(1,nrow(mat))
  return(list(Kmat=Kmat_meteo*Kmat_geno,del_theta1=Kmat_meteo_deriv*Kmat_geno,
              del_theta2=Kmat_meteo*Kmat_geno_deriv))
}


###################################################################################
# FUNCTIONS
GP_test <- function(data, train, test, betahat, Kobs_inv, Kmat, tt, vv){
  
  ind_tt <- which(data$Env == test$Env[tt] & data$variety == test$variety[tt])[1]   
  # where is the current test data in data
  Kcross <- Kmat[ind_tt, vv] %>% as.numeric()
  
  meanGP <- betahat + t(Kcross)%*%(Kobs_inv%*%(train[[tar]] - betahat*matrix(1, nrow = nrow(train), ncol = 1)))
  varGP <- Kmat[ind_tt,ind_tt] - 
    t(Kcross)%*%Kobs_inv%*%Kcross + 
    (1-t(matrix(1, nrow = nrow(train), ncol = 1))%*%Kobs_inv%*%Kcross)^2/(t(matrix(1, nrow = nrow(train), ncol = 1))%*%Kobs_inv%*%matrix(1, nrow = nrow(train), ncol = 1))
  
  pred <- meanGP
  sd <- sqrt(varGP)
  
  return(c(pred, sd))    }


select_train_ind <- function(setupp, leakage, leak_spec){
  if(setupp == '1'){
    unique_envs <- unique(data$Env)
    split_envs <- sample(unique_envs)
    train_envs <- split_envs[1:floor(length(split_envs) * train_prop)]
    test_envs <- split_envs[(floor(length(split_envs) * train_prop) + 1):length(split_envs)]}
  
  if(setupp == '2'){
    unique_vars <- unique(data$variety)
    split_vars <- sample(unique_vars)
    train_vars <- split_vars[1:floor(length(split_vars) * train_prop)]
    test_vars <- split_vars[(floor(length(split_vars) * train_prop) + 1):length(split_vars)]}
  
  if(setupp == '1' & leakage == 'no'){
    ind_train <- which(data$Env %in% train_envs)  }
  
  if(setupp == '1' & leakage == 'yes' & leak_spec == 'standard'){
    ind_train <- which(data$Env %in% train_envs | data$variety_name %in% leak_name)   }
  
  if(setupp == '1' & leakage == 'yes' & leak_spec == 'random'){
    ind_train <- which(data$Env %in% train_envs)
    # one random leak per env
    ind_leak <- NaN
    for(ttt in 1:length(test_envs)){
      vvs <- data %>% filter(Env == test_envs[ttt]) %>% pull(variety_name)
      v_rand <- sample(vvs)[1]
      ind_leak[ttt] <- which(data$Env == test_envs[ttt] & data$variety_name == v_rand)}
    ind_train <- c(ind_train, ind_leak) %>% sort()  }
  
  if(setupp == '2' & leakage == 'no'){
    ind_train <- which(data$variety_name %in% train_vars)   }
  
  if(setupp == '2' & leakage == 'yes' & leak_spec == 'standard'){
    ind_train <- which(data$variety_name %in% train_vars | data$Env %in% leak_name)   }
  
  if(setupp == '2' & leakage == 'yes' & leak_spec == 'random'){
    ind_train <- which(data$variety_name %in% train_vars)
    # one random leak per env
    ind_leak <- NaN
    for(ttt in 1:length(test_vars)){
      vvs <- data %>% filter(variety_name == test_vars[ttt]) %>% pull(Env)
      v_rand <- sample(vvs)[1]
      ind_leak[ttt] <- which(data$variety_name == test_vars[ttt] & data$Env == v_rand)}
    ind_train <- c(ind_train, ind_leak) %>% sort()  }
  
  return(ind_train)
  
}

#################################################################################
# nr_para <- 8
# registerDoParallel(cores=nr_para)
# foreach(ss = 1:nr_para)%dopar%{

log_likelihood=function(theta1,theta2,sigma,tau,k=1,train){
  if(keri == 1){
    kernel_geno <- 'ham'  
    kernel_meteo <- 'exp'   }
  if(keri == 2){
    kernel_geno <- 'spe'  
    kernel_meteo <- 'exp'   }
  if(keri == 3){
    kernel_geno <- 'ham'  
    kernel_meteo <- 'ali'   }
  if(keri == 4){
    kernel_geno <- 'spe'  
    kernel_meteo <- 'ali'   }
  
  Kmat1 <- create_kernelmat_continuous(kernel_meteo, kernel_geno,theta1,theta2,k,train)
  Kmat <- sigma*var(train[[tar]])*Kmat1
  
  if(noisy_obs==FALSE){Kobs <- Kmat }
  if(noisy_obs==TRUE){Kobs <- Kmat + tau^2*diag(nrow(train)) }
  Kobs_inv <- tryCatch(
    chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
    error = function(e) {            # Handle errors
      message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
      ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
    }
  )
  
  Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
  betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
  zn=train[[tar]]
  #logLH[i4] <- mvtnorm::dmvnorm(train[[tar]], mean = rep(betahat,nrow(train)), 
  #                              sigma = as.matrix(Kobs), log=TRUE) 
  return(nrow(train)*log(2*pi)+determinant(as.matrix(Kobs))$modulus[1]+
           t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
  
}

#mvtnorm::dmvnorm(train[[tar]], mean = rep(betahat,nrow(train)), 
#                 sigma = as.matrix(Kobs), log=TRUE) 


likelihood_gradient=function(theta1,theta2,sigma,tau,k=1,train){
  Trace=function(x){sum(diag(x))}
  if(keri == 1){
    kernel_geno <- 'ham'  
    kernel_meteo <- 'exp'   }
  if(keri == 2){
    kernel_geno <- 'spe'  
    kernel_meteo <- 'exp'   }
  if(keri == 3){
    kernel_geno <- 'ham'  
    kernel_meteo <- 'ali'   }
  if(keri == 4){
    kernel_geno <- 'spe'  
    kernel_meteo <- 'ali'   }
  
  if(kernel_geno=='spe'){
    kernel_matrices= create_kernelmat_derivative(kernel_meteo, kernel_geno,
                                                 theta1,theta2,k,train)
    Kmat1 <-kernel_matrices$Kmat
    Kmat <- sigma*var(train[[tar]])*Kmat1
    
    
    if(noisy_obs==FALSE){Kobs <- Kmat }
    if(noisy_obs==TRUE){Kobs <- Kmat + tau^2*diag(nrow(train)) }
    Kobs_inv <- tryCatch(
      chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
      error = function(e) {            # Handle errors
        message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
        ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
      }
    )
    
    Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
    betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
    zn=train[[tar]]
    
    
    Kobs_inv_deriv=Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*kernel_matrices$del_theta1)%*%Kobs_inv
    del_theta1_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    Kobs_inv_deriv=Kobs_inv%*%as.matrix(var(train[[tar]])*Kmat1)%*%Kobs_inv
    del_sigma_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    
    Kobs_inv_deriv=Kobs_inv%*%(tau*2*diag(nrow(train)))%*%Kobs_inv
    
    del_tau_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    
    grad=c(
      Trace(Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*(kernel_matrices$del_theta1)))-
        2*del_theta1_betahat%*%rep(1,nrow(train))%*%(Kobs_inv%*%(zn-betahat))-
        (zn-betahat)%*%Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*(kernel_matrices$del_theta1))%*%
        Kobs_inv%*%(zn-betahat),
      
      Trace(Kobs_inv%*%as.matrix(var(train[[tar]])*Kmat1))-
        2*del_sigma_betahat%*%rep(1,nrow(train))%*%Kobs_inv%*%(zn-betahat)-
        (zn-betahat)%*%Kobs_inv%*%as.matrix(var(train[[tar]])*Kmat1)%*%Kobs_inv%*%(zn-betahat),
      
      Trace(2*tau*Kobs_inv)-
        2*del_tau_betahat%*%rep(1,nrow(train))%*%Kobs_inv%*%(zn-betahat)-
        2*tau*(zn-betahat)%*%Kobs_inv%*%Kobs_inv%*%(zn-betahat)
      
      
      
    )
  }else{
    kernel_matrices= create_kernelmat_derivative(kernel_meteo, kernel_geno,theta1,theta2,k,train)
    Kmat1 <-kernel_matrices$Kmat
    Kmat <- sigma*var(train[[tar]])*Kmat1
    
    
    if(noisy_obs==FALSE){Kobs <- Kmat }
    if(noisy_obs==TRUE){Kobs <- Kmat + tau^2*diag(nrow(train)) }
    Kobs_inv <- tryCatch(
      chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
      error = function(e) {            # Handle errors
        message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
        ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
      }
    )
    
    Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
    betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
    zn=train[[tar]]
    
    
    Kobs_inv_deriv=Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*kernel_matrices$del_theta1)%*%Kobs_inv
    del_theta1_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    Kobs_inv_deriv=Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*kernel_matrices$del_theta2)%*%Kobs_inv
    del_theta2_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    
    Kobs_inv_deriv=Kobs_inv%*%as.matrix(var(train[[tar]])*Kmat1)%*%Kobs_inv
    del_sigma_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    
    Kobs_inv_deriv=Kobs_inv%*%(tau*2*diag(nrow(train)))%*%Kobs_inv
    
    del_tau_betahat=-sum(Kobs_inv_deriv%*%train[[tar]])/sum(Kobs_inv)+
      betahat/(sum(Kobs_inv)^2)%*%sum(Kobs_inv_deriv)
    
    
    grad=c(
      Trace(Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*(kernel_matrices$del_theta1)))-
        2*del_theta1_betahat%*%rep(1,nrow(train))%*%(Kobs_inv%*%(zn-betahat))-
        (zn-betahat)%*%Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*(kernel_matrices$del_theta1))%*%
        Kobs_inv%*%(zn-betahat),
      
      Trace(Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*(kernel_matrices$del_theta2)))-
        2*del_theta2_betahat%*%rep(1,nrow(train))%*%Kobs_inv%*%(zn-betahat)-
        (zn-betahat)%*%Kobs_inv%*%as.matrix(sigma*var(train[[tar]])*kernel_matrices$del_theta2)%*%Kobs_inv%*%(zn-betahat),
      
      Trace(Kobs_inv%*%as.matrix(var(train[[tar]])*Kmat1))-
        2*del_sigma_betahat%*%rep(1,nrow(train))%*%Kobs_inv%*%(zn-betahat)-
        (zn-betahat)%*%Kobs_inv%*%as.matrix(var(train[[tar]])*Kmat1)%*%Kobs_inv%*%(zn-betahat),
      
      Trace(2*tau*Kobs_inv)-
        2*del_tau_betahat%*%rep(1,nrow(train))%*%Kobs_inv%*%(zn-betahat)-
        2*tau*(zn-betahat)%*%Kobs_inv%*%Kobs_inv%*%(zn-betahat)
      
      
      
    )
  }
  
  return(grad)
  
}

Adam = function(f, grad_f, params_init, lr = 0.01, beta1 = 0.9, beta2 = 0.999,
                epsilon = 1e-8, max_iter = 100, tol = 1e-10,
                lb = rep(0, length(params_init)),
                ub = rep(Inf, length(params_init))) {
  
  # Initialize progress bar if not in cluster
  if (notcluster) {
    pb = progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = max_iter,
      clear = FALSE,
      width = 60
    )
  }
  
  params = params_init
  m = numeric(length(params))
  v = numeric(length(params))
  iter = 0
  
  while (iter < max_iter) {
    iter = iter + 1
    
    # Compute gradient
    grad = grad_f(params)
    
    # Update moments
    m = beta1 * m + (1 - beta1) * grad
    v = beta2 * v + (1 - beta2) * (grad^2)
    
    # Compute bias-corrected moments
    m_hat = m / (1 - beta1^iter)
    v_hat = v / (1 - beta2^iter)
    
    # Update parameters
    params_new = params - lr * m_hat / (sqrt(v_hat) + epsilon)
    
    # Handle boundary constraints
    params_new[is.na(params_new)] = params[is.na(params_new)]
    params_new[params_new <= lb] = params[params_new <= lb]
    params_new[params_new > ub] = params[params_new > ub]
    
    # Check convergence
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    params = params_new
    
    
    if (notcluster) {
      pb$tick()
    }
  }
  
  return(params)
}



for(ss in 1:2){
  
  tar <- c('yield', 'prot')[ss]
  
  hyperpp <- matrix(NA, 3, 10)
  if(keri == 1){
    kernel_geno <- 'ham'  
    kernel_meteo <- 'exp'   }
  if(keri == 2){
    kernel_geno <- 'spe'  
    kernel_meteo <- 'exp'   }
  if(keri == 3){
    kernel_geno <- 'ham'  
    kernel_meteo <- 'ali'   }
  if(keri == 4){
    kernel_geno <- 'spe'  
    kernel_meteo <- 'ali'   }
  
  tau <- mean(data[[tar]]*0.01) 
  # 1 per cent obs noise, can also be fitted as hyperparameter
  
  res <- matrix(NA,nr_outer_split,8)
  
  for(i in 1:nr_outer_split){
    
    # split train and test
    set.seed(i)
    ind_train <- select_train_ind(setupp, leakage, leak_spec)
    train <- data[ind_train,]
    test <- data[-ind_train,]
    vv <- ind_train
    
    i4 <- 1
    logLH <- NA
    for(i1 in 1:10){
      print(i1)
      for(i2 in 1:10){
        
        Kmat1 <- create_kernelmat(kernel_meteo, kernel_geno,i1,i2,train)
        
        for(i3 in 1:10){
          
          Kmat <- seq(0.1, 2.4, length.out=10)[i3]*var(data[[tar]])*Kmat1
          if(noisy_obs==FALSE){Kobs <- Kmat  }
          if(noisy_obs==TRUE){Kobs <- Kmat + tau^2*diag(nrow(train)) }
          
          Kobs_inv <- tryCatch(
            chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
            error = function(e) {            # Handle errors
              message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
              ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
            }
          )
          
          Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
          betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
          
          logLH[i4] <- mvtnorm::dmvnorm(train[[tar]], mean = rep(betahat,nrow(train)), 
                                        sigma = as.matrix(Kobs), log=TRUE) 
          i4 <- i4 + 1  }
        
      }}
    
    # if(i == 1){write.table(logLH, file = paste0("logLH", keri, tar, leakage, setupp, i, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)}
    
    mm <- which(logLH == max(logLH))
    mm1 <- ceiling(mm/100)
    mm2 <- ceiling((mm - (mm1-1)*100)/10)
    mm3 <- ifelse(mm %% 10 == 0, 10, mm %% 10)
    
    #Kmat1 <- create_kernelmat(kernel_meteo, kernel_geno,mm1,mm2,data)
    #Kmat <- seq(0.1, 2.4, length.out=10)[mm3]*var(data[[tar]])*Kmat1
    
    if(keri%%2==1){
      initial_par=c(seq(0.1, 1.9, length.out = 10)[mm1],seq(2,8,length.out=10)[mm2],
                    seq(0.1, 2.4, length.out=10)[mm3],tau)
      
      opt_par=opt_par=Adam(function(x) {
        log_likelihood(
         x[1], x[2], x[3], x[4], 1,train)
      },function(x) {
        likelihood_gradient(x[1], x[2], x[3], x[4], 1,train)
      },initial_par,lr=lr,max_iter = max_iter)
      
      Kmat1 <- create_kernelmat_continuous(kernel_meteo, kernel_geno,opt_par[1]
                                           ,opt_par[2],1,data)
      
      Kmat <- opt_par[3]*var(data[[tar]])*Kmat1
      if(noisy_obs==FALSE){Kobs <- Kmat[vv,vv]  }
      if(noisy_obs==TRUE){Kobs <- Kmat[vv,vv] + opt_par[4]^2*diag(nrow(train)) }
      
    }else{
      initial_par=c(seq(0.1, 1.9, length.out = 10)[mm1],
                    seq(0.1, 2.4, length.out=10)[mm3],tau)
      
      opt_par=opt_par=Adam(function(x) {
        log_likelihood(
          x[1], mm2, x[2], x[3], mm2,train[1:10,])
      },function(x) {
        likelihood_gradient(x[1], mm2, x[2], x[3], mm2,train[1:10,])
      },initial_par,lr=lr,max_iter = max_iter)
      
      
      Kmat1 <- create_kernelmat_continuous(kernel_meteo, kernel_geno,opt_par[1]
                                           ,mm2,mm2,data)
      
      Kmat <- opt_par[2]*var(data[[tar]])*Kmat1
      if(noisy_obs==FALSE){Kobs <- Kmat[vv,vv]  }
      if(noisy_obs==TRUE){Kobs <- Kmat[vv,vv] + opt_par[3]^2*diag(nrow(train)) }
      
      
    }
    
    Kobs_inv <- tryCatch(
      chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
      error = function(e) {            # Handle errors
        message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
        ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
      }
    )
    Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
    betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
    results <- sapply(1:nrow(test), function(tt) GP_test(data, train, test, betahat, Kobs_inv, Kmat, tt, vv))
    
    test$pred <- results[1,]
    test$sd <- sqrt(results[2,])
    # test <- test[complete.cases(test),]
    
    res[i,1] = mse(test$pred, test[[tar]])
    res[i,2] = mean(crps_norm(test[[tar]], test$pred, test$sd)) # include tau2!
    
    # global average
    mean_train <- train[[tar]] %>% mean()
    sd_train <- train[[tar]] %>% sd()
    res[i,3] = mse(test[[tar]], mean_train)
    res[i,4] = mean(abs(test[[tar]] - mean_train))
    
    # variety average
    if(leakage == 'yes' || setupp == '1'){
      pred_var <- NaN
      for(j in 1:nrow(test)){
        pred_var[j] <- train %>% 
          filter(variety_name == test$variety_name[j]) %>% 
          pull(tar) %>% mean() %>% as.numeric()}
      loca <- 1-is.nan(pred_var)
      res[i,5] = mse(test[[tar]][loca==1], pred_var[loca==1])
      res[i,6] = mean(abs(test[[tar]][loca==1] - pred_var[loca==1])) }
    
    
    # environmental average 
    if(leakage == 'yes' || setupp == '2'){
      pred_loc <- NaN
      for(j in 1:nrow(test)){
        pred_loc[j] <- train %>% ungroup() %>% 
          filter(Env == test$Env[j]) %>% 
          pull(tar) %>% mean() %>% as.numeric()}  
      loca <- 1-is.nan(pred_loc)
      res[i,7] = mse(test[[tar]][loca==1], pred_loc[loca==1])  
      res[i,8] = mean(abs(test[[tar]][loca==1] - pred_loc[loca==1]))  }
    
  }
  
  write.table(res, file = paste0('Results_hyper_', tar, keri, leakage,setupp,
                                 '.txt'), sep = "\t", row.names = FALSE, 
              quote = FALSE)
  
  
}


