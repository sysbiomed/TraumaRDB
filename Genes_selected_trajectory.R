#Among the genes previously selected by regularization with Cox model at the first 4 days 
#after injury we will represent the trajectory of the genes repeatedly selected
#Cláudia Constantino, MSc


library(ggplot2)
library(ggpubr)

#compare genes with other temporal microarrays
summary(genes_repeated_namesT0 %in% genes_repeated_namesT1)
summary(genes_repeated_namesT0 %in% genes_repeated_namesT4)
summary(genes_repeated_namesT1 %in% genes_repeated_namesT4)


genes.selected <- c(genes_repeated_namesT0,genes_repeated_namesT1, genes_repeated_namesT4)
genes.selected <- unique(genes.selected)
length(genes.selected)

genes.data.subset <- patient.data[,1]
subset <- patient.data[,genes.selected]
genes.data.subset <- cbind(genes.data.subset,subset)

#genes.data.subset has the genes select for the 797 patient observations 
summary(genes.data.subset)

genes.data.subset <- as.data.frame(genes.data.subset)


#eliminate microarray without information - Microarray ID: 23586509 and PATIENT ID: 16832704 
#correspond to index 391 in genes.data.subset dataframe
genes.data.subset <- genes.data.subset[-c(391), ] 
microarray.time <- as.data.frame(patient.data$SAMPLE_STUDYSTART_DAYS)
microarray.time <- microarray.time[-c(391),]

#get microarray time
genes.data.subset <- cbind(genes.data.subset,microarray.time)


#LINE PLOT for one of the 22 genes
ggplot(genes.data.subset, aes(x=microarray.time , 
                              y=AA808444 , 
                              group=PATIENT_ID)) +
  geom_line(color="grey45", size=0.3)+
  geom_point(color="darkslategray3", size=0.5) +
  xlab("Time [days]") + ylab("Gene 'AA808444'")


#Create 12 genes line plots in a loop
colNames <- names(genes.data.subset)[2:23]
plot_list = list()
for(i in colNames){
  n <- noquote(i)
  plt <- ggplot(genes.data.subset, aes_string(x=genes.data.subset$microarray.time , 
                                              y=genes.data.subset[,i], 
                                              group=genes.data.subset$PATIENT_ID)) +
    geom_line(color="grey45", size=0.2)+
    geom_point(color="darkslategray3", size=0.1) +
    xlab("Time [days]") + ylab(n)
  print(plt)
  plot_list[[i]] = plt
  Sys.sleep(2)
}

figure1 <- ggarrange(plotlist = plot_list, ncol = 6, nrow = 4)
figure1


#########