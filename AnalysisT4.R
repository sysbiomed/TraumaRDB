#Analysis of gene data at the fourth day post injury (111 patients included)
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



#select in data t0 for a simple cox analysis
data_t4 <- patient.data[(patient.data$SAMPLE_STUDYSTART_DAYS == "4"),]
data_t4$SAMPLE_STUDYSTART_DAYS <- NULL

#create vector with 1 and 0 corresponding to event
colnames(data_t4)[3] <- "EVENT"
data_t4$EVENT <- ifelse(data_t4$EVENT > 1 , 0, data_t4$EVENT)
data_t4$EVENT[is.na(data_t4$EVENT)] <- 1 #substitue NA's for 0


ydata_t4 <- cbind(data_t4$DISCHARGE_DAY_SNC_INJ, data_t4$EVENT)
colnames(ydata_t4) <- c("time", "status")

genes_t4 <- as.data.frame(data_t4[c(5:54769)])
genes_t4 <- as.matrix(genes_t4)



# =====================================================================================================
# See if possible to predict time of discharge using regularization cox for T0
# =====================================================================================================
#TESTS WITH DATA SPLIT INTO TEST/TRAIN

#split 70% for train and 30% for test 
index <- sample(nrow(genes_t4), 78)
datatreino <- as.matrix(genes_t4[index, ])
datateste <- as.matrix(genes_t4[-index, ])
vectortreino <- ydata_t4[index, ]
vectorteste <- ydata_t4[-index, ]


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
                   plot.title = 'Train set at T4', legend.outside = FALSE)


#TEST
vectorteste <- as.data.frame(vectorteste)
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Test set at T4', legend.outside = FALSE)


# =====================================================================================================
# SELECT GENES REPEATEDLY ASSOCIATED WITH THE EVENT
# =====================================================================================================
#LEAVE ONE OUT METHOD 
ydata_t4 <- as.matrix(ydata_t4)
#loo <- NULL
loo_names <- NULL
degree = 1:111 #maximum = 111

for (d in degree){
  datan <- genes_t4[-d,]
  response.vectorn <- ydata_t4[-d,]
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

loo_names_111 <- loo_names #2675 genes


loo_names <- loo_names_111
#To count the number of repeated genes in the 100 models
#with method Leave One Out 
loo_namest <- as.vector(t(loo_names))
loo_selectedgenes <- table(loo_namest)
loo_selectedgenes <- as.data.frame(loo_selectedgenes)
totalgenes_selected <- as.data.frame(loo_selectedgenes$Freq)
totalgenes_selected <- as.vector(t(totalgenes_selected)) #NR OF GENES ALWAYS SELECTED
table(totalgenes_selected)


repeated.genesT4 <- filter(loo_selectedgenes, Freq > 100) #90% 
allrepeated.genesT4 <- filter(loo_selectedgenes, Freq == 111) #100%



#if we use just this genes repeated for model prediction:
genes_repeated_namesT4 <- as.vector(allrepeated.genesT4$loo_namest)
genesT4_repeated <- as.matrix(genes_t4[,genes_repeated_namesT4])

#split 70% for train and 30% for test 
index <- sample(nrow(genesT4_repeated), 78)
datatreino <- as.matrix(genesT4_repeated[index, ])
datateste <- as.matrix(genesT4_repeated[-index, ])
vectortreino <- ydata_t4[index, ]
vectorteste <- ydata_t4[-index, ]


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
                   plot.title = 'Train set at T4 - model with repeatedly selected genes in LOOCV',
                   legend.outside = FALSE)


#TEST
vectorteste <- as.data.frame(vectorteste)
separate2GroupsCox(as.vector(elastic.cox), 
                   datateste[, names(elastic.cox)], 
                   vectorteste, 
                   plot.title = 'Test set at T4 - model with repeatedly selected genes in LOOCV',
                   legend.outside = FALSE)


#END
