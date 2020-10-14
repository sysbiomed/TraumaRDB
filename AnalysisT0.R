#Analysis at day 0 after injury, as time independent variables
#Cláudia Constantino, MSc 


if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("glmSparseNet")

library(dplyr)
library(ggplot2)
library(survival)
library(loose.rock)
library(futile.logger)
#
library(glmSparseNet)
library(glmnet)
#
library(tidyverse)
library(readxl)



#select in data t0 for a survival analysis
data_t0 <- patient.data[(patient.data$SAMPLE_STUDYSTART_DAYS == "0"),]
data_t0$SAMPLE_STUDYSTART_DAYS <- NULL

colnames(data_t0)[3] <- "EVENT"
data_t0$EVENT <- ifelse(data_t0$EVENT > 1 , 0, data_t0$EVENT)
data_t0$EVENT[is.na(data_t0$EVENT)] <- 1 #substitue NA's for 0

ydata_t0 <- cbind(data_t0$DISCHARGE_DAY_SNC_INJ, data_t0$EVENT)
colnames(ydata_t0) <- c("time", "status")


genes_t0 <- as.matrix(data_t0[c(5:54679)])
ydata_t0 <- as.matrix(ydata_t0)



# =====================================================================================================
# See if possible to predict time of discharge using regularization cox for T0
# =====================================================================================================
#TESTS WITH DATA SPLIT INTO TEST/TRAIN

#split 70% for train and 30% for test 
index <- sample(nrow(genes_t0), 118)
datatreino <- as.matrix(genes_t0[index, ])
datateste <- as.matrix(genes_t0[-index, ])
vectortreino <- ydata_t0[index, ]
vectorteste <- ydata_t0[-index, ]

vectortreino <- as.matrix(vectortreino)
datatreino <- as.matrix(datatreino)
cvfit = cv.glmnet(datatreino, vectortreino, family = "cox")


fit <- glmnet(datatreino, vectortreino, 
              alpha = 0.8,
              lambda = cvfit$lambda.min,
              family = "cox")

elastic.cox <- coef(fit, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select just variables with coef different from zero
length(elastic.cox)
names(elastic.cox)


#to compute survival curve
#TREINO
vectortreino <- as.data.frame(vectortreino)
separate2GroupsCox(as.vector(elastic.cox), 
                   datatreino[, names(elastic.cox)], 
                   vectortreino, 
                   plot.title = 'Train set at T0', legend.outside = FALSE)


#TEST
vectorteste <- as.data.frame(vectorteste)
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Test set at T0', legend.outside = FALSE)





# =====================================================================================================
# SELECT GENES REPEATEDLY ASSOCIATED WITH THE EVENT
# =====================================================================================================
#LEAVE ONE OUT METHOD 
ydata_t0 <- as.matrix(ydata_t0)
#loo <- NULL
loo_names <- NULL
degree = 91:168

for (d in degree){
  datan <- genes_t0[-d,]
  response.vectorn <- ydata_t0[-d,]
  cvfit = cv.glmnet(datan, response.vectorn, family = "cox")
  
  model_loo <- glmnet(datan,
                      response.vectorn,
                      alpha = 0.8,
                      lambda = cvfit$lambda.min,
                      family = "cox")
  
  loo_elastic <- coef(model_loo, s = 'lambda.min')[,1] %>% {.[.!=0]}
  length(loo_elastic)
  loo_elasticnames <- as.data.frame(names(loo_elastic)) #ou as.vector
  loo_names <- rbind(loo_names, loo_elasticnames)
  
  Sys.sleep(0.1)
  print(d)
}

loo_names_168 <- loo_names #2919 genes



#To count the number of repeated genes in the 168 models
#with method Leave One Out 
loo_namest <- as.vector(t(loo_names))
loo_selectedgenes <- table(loo_namest)
loo_selectedgenes <- as.data.frame(loo_selectedgenes)
totalgenes_selected <- as.data.frame(loo_selectedgenes$Freq)
totalgenes_selected <- as.vector(t(totalgenes_selected)) #NR OF GENES ALWAYS SELECTED
table(totalgenes_selected)


repeated.genesT0 <- filter(loo_selectedgenes, Freq > 150) #90%
allrepeated.genesT0 <- filter(loo_selectedgenes, Freq == 168) #100%


#if we use just this genes repeated for model prediction:
genes_repeated_namesT0 <- as.vector(allrepeated.genesT0$loo_namest)
genesT0_repeated <- as.matrix(genes_t0[,genes_repeated_namesT0])

#split 70% for train and 30% for test 
index <- sample(nrow(genesT0_repeated), 118)
datatreino <- as.matrix(genesT0_repeated[index, ])
datateste <- as.matrix(genesT0_repeated[-index, ])
vectortreino <- ydata_t0[index, ]
vectorteste <- ydata_t0[-index, ]


vectortreino <- as.matrix(vectortreino)
cvfit = cv.glmnet(datatreino, vectortreino, family = "cox")

fit <- glmnet(datatreino, vectortreino, 
              alpha = 0,
              lambda = cvfit$lambda.min,
              family = "cox")

elastic.cox <- coef(fit, s = 'lambda.1se')[,1] %>% {.[.!=0]} #to select just variables with coef different from zero
length(elastic.cox)
names(elastic.cox)


#to compute survival curve
#TREINO
vectortreino <- as.data.frame(vectortreino)
separate2GroupsCox(as.vector(elastic.cox), 
                   datatreino[, names(elastic.cox)], 
                   vectortreino, 
                   plot.title = 'Train set at T0 - model with repeatedly selected genes in LOOCV', legend.outside = FALSE)


#TEST
vectorteste <- as.data.frame(vectorteste)
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Test set at T0 - model with repeatedly selected genes in LOOCV', legend.outside = FALSE)


#END



