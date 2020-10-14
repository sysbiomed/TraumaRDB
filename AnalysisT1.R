#Analysis of gene data at the day 1 post injury (162 patients included)
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
library(tidyverse)
library(readxl)



#select in data t0 for a survival analysis
data_t1 <- patient.data[(patient.data$SAMPLE_STUDYSTART_DAYS == "1"),]
data_t1$SAMPLE_STUDYSTART_DAYS <- NULL

#create vector with 1 and 0 corresponding to event
colnames(data_t1)[3] <- "EVENT"
data_t1$EVENT <- ifelse(data_t1$EVENT > 1 , 0, data_t1$EVENT)
data_t1$EVENT[is.na(data_t1$EVENT)] <- 1 #substitue NA's for 0


ydata_t1 <- cbind(data_t1$DISCHARGE_DAY_SNC_INJ, data_t1$EVENT)
colnames(ydata_t1) <- c("time", "status")

genes_t1 <- as.data.frame(data_t1[c(5:54679)])
genes_t1 <- as.matrix(genes_t1)




# =====================================================================================================
# See if possible to predict time of discharge using regularization cox for T0
# =====================================================================================================
#TESTS WITH DATA SPLIT INTO TEST/TRAIN

#split 70% for train and 30% for test 
index <- sample(nrow(genes_t1), 113)
datatreino <- as.matrix(genes_t1[index, ])
datateste <- as.matrix(genes_t1[-index, ])
vectortreino <- ydata_t1[index, ]
vectorteste <- ydata_t1[-index, ]


vectortreino <- as.matrix(vectortreino)
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
                   plot.title = 'Train set at T1', legend.outside = FALSE)


#TEST
vectorteste <- as.data.frame(vectorteste)
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Test set at T1', legend.outside = FALSE)




# =====================================================================================================
# SELECT GENES REPEATEDLY ASSOCIATED WITH THE EVENT
# =====================================================================================================
#LEAVE ONE OUT METHOD 
ydata <- as.matrix(ydata)
#loo <- NULL
loo_names <- NULL
degree = 91:162 #maximum = 162

for (d in degree){
  datan <- genes_t1[-d,]
  response.vectorn <- ydata[-d,]
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

loo_names_162 <- loo_names #2675 genes


#To count the number of repeated genes in the 100 models
#with method Leave One Out 
loo_namest <- as.vector(t(loo_names))
loo_selectedgenes <- table(loo_namest)
loo_selectedgenes <- as.data.frame(loo_selectedgenes)
totalgenes_selected <- as.data.frame(loo_selectedgenes$Freq)
totalgenes_selected <- as.vector(t(totalgenes_selected)) #NR OF GENES ALWAYS SELECTED
table(totalgenes_selected)


repeated.genesT1 <- filter(loo_selectedgenes, Freq > 145) #90% 
allrepeated.genesT1 <- filter(loo_selectedgenes, Freq == 162) #100%
allrepeated.genesT1 <- allrepeated.genesT1[-9,] #apagar esta linha


#if we use just this genes repeated for model prediction:
genes_repeated_namesT1 <- as.vector(allrepeated.genesT1$loo_namest)
genesT1_repeated <- as.matrix(genes_t1[,genes_repeated_namesT1])

#split 70% for train and 30% for test 
index <- sample(nrow(genesT1_repeated), 113)
datatreino <- as.matrix(genesT1_repeated[index, ])
datateste <- as.matrix(genesT1_repeated[-index, ])
vectortreino <- ydata[index, ]
vectorteste <- ydata[-index, ]


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
                   plot.title = 'Train set at T1 - model with repeatedly selected genes in LOOCV',
                   legend.outside = FALSE)


#TEST
vectorteste <- as.data.frame(vectorteste)
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Test set at T1 - model with repeatedly selected genes in LOOCV',
                   legend.outside = FALSE)


#END









