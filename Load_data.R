#Import large GlueGrant gene expression data
#Cláudia Constantino, MSc 


library(readxl)
library(dplyr)


#Import clinical variables that will be used in the study
clinical.data <- read_excel(
  "./data/clinical_data/TRDB_TRAUMA-PATIENT_CLIN_DATA_RPT_20090126_195655.xls")
clinical.data <- subset(clinical.data, select = c(PATIENT_ID, DISCHARGE_DAY_SNC_INJ, DEATH_DAY_SNC_INJ))
head(clinical.data)


#Import the correspondence between patient and microarray
patient.microarray <- read_excel(
  "./data/clinical_data/TRDB_TRAUMA-PT_DEMO_MICRO_SVRTY_RPT_20090824_105202.xls")
patient.microarray <- subset(patient.microarray, select = c(PATIENT_ID, MICROARRAY_ID,
                                                            SAMPLE_STUDYSTART_DAYS))
head(patient.microarray)

#Inspect the microarray data
microarray <- read_excel("./data/normarray/microarrayIDALL.168.batch1/microarrayIDALL.168.batch1.expression.xls")
microarray <- as.data.frame(microarray)
rownames(microarray) <- microarray$`probe set`
microarray <- microarray[, (names(microarray) != "probe set")]
head(microarray[, 1:10])



# Prealocate data.frame
microarray.data <- as.data.frame(matrix(NA, nrow(patient.microarray), nrow(microarray)+1))
colnames(microarray.data) <- c("MICROARRAY_ID", rownames(microarray))
microarray.data$MICROARRAY_ID <- patient.microarray$MICROARRAY_ID



#import the remaining folders with gene expression information
excel.data1 <- read_excel("./data/normarray/microarrayIDALL.168.batch1/microarrayIDALL.168.batch1.expression.xls")
excel.data1 <- as.data.frame(excel.data1)
colnames(excel.data1)[1] <- "probe.set"
excel.data2 <- read.csv("./data/normarray/microarrayIDALL.168.batch2/microarrayIDALL.168.batch2.expression.csv")
excel.data3 <- read.csv("./data/normarray/microarrayIDALL.168.batch3/microarrayIDALL.168.batch3.expression.csv")
excel.data4 <- read.csv("./data/normarray/microarrayIDALL.168.batch4/microarrayIDALL.168.batch4.expression.csv")

my_list <- list(excel.data1, excel.data2, excel.data3, excel.data4)

for(i in 1:length(my_list)){
  excel.data <- as.data.frame(my_list[i])
  
  row.names(excel.data) <- excel.data$probe.set
  excel.data <- excel.data[, names(excel.data) != "probe.set"]
  
  for(i in 1:ncol(excel.data)){
    col.name <- colnames(excel.data)[i]
    microarray.number <- as.numeric(sub("t(\\d+)", "\\1", col.name))
    
    microarray.data[microarray.data$MICROARRAY_ID == microarray.number, 
                    -1] <- excel.data[, i]
  }
}
microarray.data[1:10, 1:10]




#Substitute the probe ID into the accession number (genbank reference number)
probeID.to.accession <- read_excel("./ProbeID_to_Accession.xlsx")
probeID <- data.frame("ProbeID" = colnames(microarray.data)[-1])
accession <- inner_join(probeID, probeID.to.accession, by = "ProbeID")


# The probes that start with AFFX don't have a accession number. We'll use the ProbeID
accession[is.na(accession$Accession), "Accession"] <- accession[is.na(accession$Accession), "ProbeID"] 
probes.with.accession <- accession[!is.na(accession$Accession), 1]
colnames(microarray.data) <- c("MICROARRAY_ID", accession[!is.na(accession$Accession), 2])
microarray.data[1:10, 1:10]


#add data to patients id 
patient.with.microarray <- inner_join(clinical.data, patient.microarray, by = "PATIENT_ID", all.x = TRUE)
#cofirm that the order of the microarrays are equal in both dataframes
identical(patient.with.microarray[['MICROARRAY_ID']],microarray.data[['MICROARRAY_ID']])

patient.data <- cbind(patient.with.microarray, microarray.data[,2:54676])
View(patient.with.microarray)




