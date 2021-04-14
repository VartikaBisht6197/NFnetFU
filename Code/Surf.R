#set working directory
setwd("/Users/Vartika_Bisht/Individual_Project")

library(Jmisc)

#dowload all files from http://www.mathstat.dal.ca/tkenney/Rpackages/
sourceAll("/Users/Vartika_Bisht/Individual_Project/SuRF/R")

# Dataset 1 SuRF
Data_Set_1 <- read.xlsx("NFnetFU_Dataset1_wih_labels.xlsx", sheetIndex = 1)
df_data1 <- as.data.frame(Data_Set_1[,3:49])
data1 <- data.matrix(df_data1)
X <- data1
Y <- as.numeric(factor(Data_Set_1$Class))
fitting<-SURF(X,Y,display.progress=TRUE,family=stats::poisson(link="log"))

# Dataset 2 SuRF
Data_Set_1 <- read.csv("/Users/Vartika_Bisht/NFnetFU/Data/NFnetFU_Dataset2_OTU_abundance(top 100).csv")
df_data1 <- as.data.frame(Data_Set_1)
data1 <- data.matrix(df_data1[,2:101])
X <- data1
meta <- read.csv("/Users/Vartika_Bisht/NFnetFU/Data/NFnetFU_Dataset2_metadata.csv")
Y <- as.numeric(factor(meta$dx))
fitting2<-SURF(X,Y,display.progress=TRUE,family=stats::poisson(link="log"))


# Dataset 3 SuRF
Data_Set_1 <- read.csv("/Users/Vartika_Bisht/NFnetFU/Data/NFnetFU_Dataset3_OTU_abundance(top 100).csv")
df_data1 <- as.data.frame(Data_Set_1)
data1 <- data.matrix(df_data1[,2:101])
X <- data1
meta <- read.csv("/Users/Vartika_Bisht/NFnetFU/Data/NFnetFU_Dataset3_metadat.csv")
Y <- as.numeric(factor(meta$chem_administration[1:422]))
fitting3<-SURF(X,Y,display.progress=TRUE,family=stats::poisson(link="log"))

