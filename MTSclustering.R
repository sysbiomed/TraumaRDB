#Clustering of multivariate time series (MTS clustering) with "dtwclust" package
#Cláudia Constantino, MSc 

install.packages("dtwclust")
library(dtwclust)
library (cluster)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)



# ====================================================================================
#CLUSTERING OF MULTIVARIATE TIME SERIES with completed patients info (n = 33)
# ====================================================================================

#to get the necessary arranged input dataframe to enter in the method: 
# "Multivariate series should be provided as a list of matrices where time spans 
#  the rows and the variables span the columns of each matrix."
genes22info <- bind_rows(ls_imput_complete)

genes22info <- as.data.frame(genes22info)
genes22info <- genes22info %>% arrange(patient)


#create a list
list_genes22info <- split(genes22info , f = genes22info$patient )
list_genes22info <- lapply(list_genes22info, function(x) { x["patient"] <- NULL; x })

#Transpose each data frame in the list (days in rows and genes in columns)
Tlist_genes22info <- lapply(list_genes22info,function(x){ t(x) })



#to get the better p-value for different k
#get the time and event information just for the 33 patients included in this step
cox.data.genesclust <- as.data.frame(ls_imput_complete[["AA808444"]]$patient)

df.time.event <- cbind(data_t0$PATIENT_ID , ydata_t0)
df.time.event <- as.data.frame(df.time.event)

ids <- as.vector(cox.data.genesclust[,1])
ids

df.time.event.33 <- df.time.event[df.time.event$V1 %in% ids,]

cox.data.genesclust <- cbind(cox.data.genesclust,
                             df.time.event.33)
cox.data.genesclust <- as.data.frame(cox.data.genesclust)
colnames(cox.data.genesclust) <- c("patient" , "patient2", "time", "status")


#MTS clustering for a range of K and Log Rank test
#using function tsclust from "dtwclust" package 
for (i in 5:10){
  mvclustering <- tsclust(Tlist_genes22info, k = i, distance = "GAK", seed = 390,
                          args = tsclust_args(dist = list(sigma = 100)))
  
  cox.data.genesclust[,5] <- mvclustering@cluster
  colnames(cox.data.genesclust)[5] <- c("MTScluster")
  
  
  #Log-Rank test comparing survival curves
  surv_diff <- survdiff(Surv(cox.data.genesclust$time, cox.data.genesclust$status) ~ 
                          MTScluster, 
                        data = cox.data.genesclust)
  print(surv_diff)
  
}


#with the k choosen in the previous step
mvclustering <- tsclust(Tlist_genes22info, k = 10L, distance = "GAK", seed = 390,
                        args = tsclust_args(dist = list(sigma = 100)))
# Note how the variables of each series are appended one after the other in the plot
plot(mvclustering)


#univariate cox analysis for the MTS clustering with K=8 
cox.data.genesclust[,5] <- mvclustering@cluster
colnames(cox.data.genesclust)[5] <- c("MTScluster")
summary(cox.data.genesclust)
cox.data.genesclust$MTScluster <- as.character(cox.data.genesclust$MTScluster)

cox.genesclust <- coxph(Surv(cox.data.genesclust$time, cox.data.genesclust$status) 
                        ~ MTScluster , data = cox.data.genesclust)
summary(cox.genesclust)


#Kaplan Meier Analysis
Kaplan.surv <- survfit(Surv(cox.data.genesclust$time, cox.data.genesclust$status) ~ 
                         MTScluster, 
                       data = cox.data.genesclust)
print(Kaplan.surv)


ggsurvplot(Kaplan.surv,
           pval = TRUE, conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

#Log-Rank test comparing survival curves
surv_diff <- survdiff(Surv(cox.data.genesclust$time, cox.data.genesclust$status) ~ 
                        MTScluster, 
                      data = cox.data.genesclust)
print(surv_diff)







# ====================================================================================
#CLUSTERING OF MULTIVARIATE TIME SERIES with completed patients until T7 info (n = 141)
# ====================================================================================

#to get the necessary arranged input dataframe to enter in the method: 
# "Multivariate series should be provided as a list of matrices where time spans 
#  the rows and the variables span the columns of each matrix."
genes22_141ids <- bind_rows(ls_imput_T7)

genes22_141ids <- as.data.frame(genes22_141ids)
genes22_141ids <- genes22_141ids %>% arrange(patient)


#create a list
list_genes22_T7 <- split(genes22_141ids , f = genes22_141ids$patient )
list_genes22_T7 <- lapply(list_genes22_T7, function(x) { x["patient"] <- NULL; x })

#Transpose each data frame in the list (days in rows and genes in columns)
Tlist_genes22_T7 <- lapply(list_genes22_T7,function(x){ t(x) })


#to get the better p-value for different k
#get the time and event information just for the 33 patients included in this step
cox.data.genesclust.T7 <- as.data.frame(ls_imput_T7[["AA808444"]]$patient) 

df.time.event.T7 <- cbind(data_t0$PATIENT_ID , ydata_t0)
df.time.event.T7 <- as.data.frame(df.time.event.T7)

ids.T7 <- as.vector(cox.data.genesclust.T7[,1])
ids.T7

df.time.event.T7 <- df.time.event.T7[df.time.event.T7$V1 %in% ids.T7,]

cox.data.genesclust.T7 <- cbind(cox.data.genesclust.T7,
                                df.time.event.T7)
cox.data.genesclust.T7 <- as.data.frame(cox.data.genesclust.T7)
colnames(cox.data.genesclust.T7) <- c("patient" , "patient2", "time", "status")


#MTS clustering with for a range of K and Log Rank test
#using function tsclust from "dtwclust" package 
for (i in 10:20){
  mvclustering_T7 <- tsclust(Tlist_genes22_T7, k = i, distance = "GAK", seed = 390,
                             args = tsclust_args(dist = list(sigma = 100)))
  
  cox.data.genesclust.T7[,5] <- mvclustering_T7@cluster
  colnames(cox.data.genesclust.T7)[5] <- c("MTScluster")
  
  cox.data.genesclust.T7$MTScluster <- as.character(cox.data.genesclust.T7$MTScluster)
  
  cox.genesclust.T7 <- coxph(Surv(cox.data.genesclust.T7$time, cox.data.genesclust.T7$status) 
                             ~ MTScluster , data = cox.data.genesclust.T7)
  summary(cox.genesclust.T7)
  
  #Log-Rank test comparing survival curves
  surv_diff <- survdiff(Surv(cox.data.genesclust.T7$time, cox.data.genesclust.T7$status) ~ 
                          MTScluster, 
                        data = cox.data.genesclust.T7)
  print(surv_diff)
  
}


#with the k chosen in the previous step
mvclustering_T7 <- tsclust(Tlist_genes22_T7, k = 15L, distance = "GAK", seed = 390,
                           args = tsclust_args(dist = list(sigma = 100)))
# Note how the variables of each series are appended one after the other in the plot
plot(mvclustering_T7)



#univariate cox analysis for the MTS clustering with K chosen 
cox.data.genesclust.T7[,5] <- mvclustering_T7@cluster
colnames(cox.data.genesclust.T7)[5] <- c("MTScluster")
summary(cox.data.genesclust.T7)
cox.data.genesclust.T7$MTScluster <- as.character(cox.data.genesclust.T7$MTScluster)

cox.genesclust.T7 <- coxph(Surv(cox.data.genesclust.T7$time, cox.data.genesclust.T7$status) 
                           ~ MTScluster , data = cox.data.genesclust.T7)
summary(cox.genesclust.T7)



#Kaplan Meier Analysis
Kaplan.surv <- survfit(Surv(cox.data.genesclust.T7$time, cox.data.genesclust.T7$status) ~ 
                         MTScluster, 
                       data = cox.data.genesclust.T7)
print(Kaplan.surv)


ggsurvplot(Kaplan.surv,
           pval = TRUE, conf.int = FALSE,
           risk.table = FALSE, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)


#Log-Rank test comparing survival curves
surv_diff <- survdiff(Surv(cox.data.genesclust.T7$time, cox.data.genesclust.T7$status) ~ 
                        MTScluster, 
                      data = cox.data.genesclust.T7)
print(surv_diff)


#####

