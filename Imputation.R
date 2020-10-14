#Automated code used for the new imputation method (junction of LOCF and NOCB)
#Cláudia Constantino, MSc 

install.packages("xts", repos="http://cloud.r-project.org")
install.packages("imputeTS")

library(dplyr)
library(tidyr)
library(zoo)
library(xts)
library(lubridate)
library(imputeTS)
library(ggplot2)
library(ggpubr)
library(Rcpp)



#Log transformation in the subset with 22 genes
Log.aux <- log1p(genes.data.subset[,2:23])
Log.genes.subset <- cbind(genes.data.subset$PATIENT_ID, 
                          genes.data.subset$microarray.time,
                          Log.aux)

colnames(Log.genes.subset)[1:2] <- c("patient", "time")

# ====================================================================================
# Carry Forward and Carry Backward Imputation
# ====================================================================================

NamesList <- colnames(Log.genes.subset)
NamesList <- NamesList[-(1:2)]
geneX <- Log.genes.subset[,1:2]
ls_imput<-list()


for (j in 1:length(NamesList)){
  geneX <- cbind(Log.genes.subset[,1:2], Log.genes.subset[NamesList[j]])
  
  #columns divided by time 
  wide_geneX <- geneX %>% spread(time, NamesList[j])
  colnames(wide_geneX)[2:23] <- paste("Time", colnames(wide_geneX[,c(2:23)]), sep = "")
  
  #Imputation: Carry near observation (controlled)
  
  #Carry the t+1 observation to t observation
  new_imput_geneX = within(wide_geneX, { 
    Time1 = ifelse(is.na(wide_geneX$Time1), wide_geneX$Time2, wide_geneX$Time1) 
    Time4 = ifelse(is.na(wide_geneX$Time4), wide_geneX$Time5, wide_geneX$Time4)
    Time7 = ifelse(is.na(wide_geneX$Time7), wide_geneX$Time8, wide_geneX$Time7)
    Time14 = ifelse(is.na(wide_geneX$Time14), wide_geneX$Time15, wide_geneX$Time14)
    Time21 = ifelse(is.na(wide_geneX$Time21), wide_geneX$Time22, wide_geneX$Time21)
    Time28 = ifelse(is.na(wide_geneX$Time28), wide_geneX$Time29, wide_geneX$Time28)
  } )
  
  #Carry the t-1 observation to t observation
  new_imput_geneX = within(new_imput_geneX, { 
    Time4 = ifelse(is.na(new_imput_geneX$Time4), new_imput_geneX$Time3, new_imput_geneX$Time4)
    Time7 = ifelse(is.na(new_imput_geneX$Time7), new_imput_geneX$Time6, new_imput_geneX$Time7)
    Time14 = ifelse(is.na(new_imput_geneX$Time14), new_imput_geneX$Time13, new_imput_geneX$Time14)
    Time21 = ifelse(is.na(new_imput_geneX$Time21), new_imput_geneX$Time20, new_imput_geneX$Time21)
    Time28 = ifelse(is.na(new_imput_geneX$Time28), new_imput_geneX$Time27, new_imput_geneX$Time28)
  } )
  
  
  #Carry the t+2 observation to t observation
  new_imput_geneX = within(new_imput_geneX, { 
    Time7 = ifelse(is.na(new_imput_geneX$Time7), new_imput_geneX$Time9, new_imput_geneX$Time7)
    Time14 = ifelse(is.na(new_imput_geneX$Time14), new_imput_geneX$Time16, new_imput_geneX$Time14)
    Time21 = ifelse(is.na(new_imput_geneX$Time21), new_imput_geneX$Time23, new_imput_geneX$Time21)
  } )
  
  
  #delete days that will not be used 
  new_imput_geneX[ ,c('Time2', 'Time3', 'Time5', 'Time6', 'Time8', 'Time9', 'Time11', 'Time13',
                      'Time15', 'Time16', 'Time20', 'Time22', 'Time23', 'Time27',
                      'Time29')] <- list(NULL)
  
  
  #And now, to complete the remaining NA's, linear interpolation
  colnames(new_imput_geneX) <- c("patient", "0", "1", "4", "7", "14", "21", "28")
  long_geneX_imput1 <- gather(new_imput_geneX, time, geneX,"0":"28")
  long_geneX_imput1 <- transform(long_geneX_imput1, time = as.numeric(time))
  long_geneX_imput1 <- long_geneX_imput1[order(long_geneX_imput1$patient,
                                               long_geneX_imput1$time),]
  #to do linear interpolation for each patient and not the entire column
  wide_aux <- long_geneX_imput1 %>% spread(patient, geneX)
  
  
  wide_geneX_imput2 <- as.data.frame(wide_aux$time)
  
  #vector with gene expression gene data to enter for the imputation
  for (i in names(wide_aux)) {
    v_aux <- as.vector(wide_aux[[i]])
    index <- max(which(!is.na(v_aux))) +1
    if (index <= 7) {
      v_aux <- v_aux[-c(index:7)]
      imput2_aux <- na_interpolation(v_aux, option = "linear")
      imput2_aux[c(index:7)] <- NA
    }
    else{
      imput2_aux <- na_interpolation(v_aux, option = "linear")
    }
    wide_geneX_imput2[,i] <- as.data.frame(imput2_aux)
  }
  wide_geneX_imput2[,1] <- NULL
  
  long_geneX_imput2 <- gather(wide_geneX_imput2, patient, geneX ,"37134":"34912227")
  new_imput2_geneX <- long_geneX_imput2 %>% spread(time, geneX) #wide
  new_imput2_geneX <- new_imput2_geneX[ order(match(new_imput2_geneX$patient, 
                                                    new_imput_geneX$patient)), ]
  
  ls_imput[[NamesList[j]]] <-new_imput2_geneX #list which contains the dataframes with imputation for each gene
  
}


#how many patients have no missing values in each time point after this imputation
ls_imput_complete <- lapply(ls_imput, na.omit)

#how many patients have no missing values until T7 after this imputation
ls_imput_T7 <- lapply(ls_imput, function(x) { x[6:8] <- NULL; x })
ls_imput_T7 <- lapply(ls_imput_T7, na.omit)








