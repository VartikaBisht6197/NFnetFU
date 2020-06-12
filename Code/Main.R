setwd("/Users/Vartika_Bisht/Individual_Project")
source("Penalty_Function.R")
source("Incorporate_Groups.R")
source("Required_Libraries.R")
source("Borrowed_Functions.R")
source("MicrobiomeAnalyst.R")

#Input csv
Data_Set_1 <- read.csv("adenoma.csv")[,2:101]
meta_data <- read.csv("meta date.csv")

#Change the file in the format required by Module 1
label_dat <- as.numeric(factor(meta_data$dx))
data1 <- data.matrix(Data_Set_1)

##-----Data Driven Module-----

## Input : Numeric Labels(label_dat) and Microbiome Abundance Data(data1)
source("Module_1.R")
## Output : Rule Based Matrix (rules_int) , Scaled Rule Based Matrix (scaled_rules_int) and Labels (label_dat)

## Input : Scaled Rule Based Matrix (scaled_rules_int)
source("Module_2.R")
## Output : Rule Based matrix with Colinearity Handled (new_data1) and PCA Loadings used to combine groups (PCA_loadings) 

## Input : Rule Based matrix with Colinearity Handled (new_data1) and PCA Loadings used to combine groups (PCA_loadings) 
source("Module_3.R")
## Output : Feature Parameters (feature_parameters)


##-----Network Fusion and Visualisation-----

#Diseases to look for in TSEA
disease <- c("Colorectal","Crohn")

#List of Microboes from selected features (OTU)
#OTU to Microbes
OTU_file <- read.table("final.csv", header = 1)
OTU_index <- which(OTU_file$OTU %in% rownames(feature_parameters))
selected_OTU <- OTU_file[OTU_index,]
feature_inorder <- selected_OTU$OTU
write.csv(selected_OTU,"OTU Microbes Selected Table.csv")

#Valid Microbe Names
OTU_network <- c()
taxa <- strsplit(as.character(selected_OTU$Taxonomy),";")
for(i in 1:length(taxa)){
  taxa_name <- taxa[[i]][length(taxa[[i]])]
  OTU_network <- c(OTU_network, substr( taxa_name , 1 , nchar(taxa_name)-5) )
}


#List of Microboes from selected features (Microbes)
Microbes_name <- substring(colnames(rules_int),4)
OTU_network <- c()
for(i in Microbes_name){
  n <- strsplit(i,split='.', fixed=TRUE)[[1]]
  if((length(n)>1)&&(n[2] == "unidentified")){
    OTU_network <- c(OTU_network,sprintf("%s.%s",n[1],n[2]))
  }else{
    OTU_network <- c(OTU_network,n[1]) 
  }
}
feature_inorder <- colnames(rules_int)

Microbes <- unique(OTU_network)


## Input : List of Microbes
source("Module_4(TSEA Network).R")
## Output : Network and Network Legends with Node size (Legends)


Name_Change <- as.data.frame(OTU_network)
rownames(Name_Change) <- feature_inorder
write.csv(Name_Change,"Features to Microbes for TSEA.csv")


## Input : Infusing Data Driven Information
Cluster_Parameters <- c()
Cluster_OTU_name <- c()
for (i in Network_Info[,"Microbe Names"]) {
  index <- which(OTU_network %in% i)
  OTU <- as.character(feature_inorder[index])
  Cluster_OTU_name <- append(Cluster_OTU_name,list(OTU))
  if(length(OTU) > 1){
    OTUs_val <- c()
    for (j in OTU) {
      OTUs_val <- c( OTUs_val , abs(feature_parameters[j,]) )
    }
    CP <- (sum(OTUs_val)/length(OTUs_val))[1]
  } else {
    CP <- abs(feature_parameters[OTU,])
  }
  Cluster_Parameters <- c( Cluster_Parameters , CP )
}
Data_Bio_Driven <- cbind(Network_Info,Cluster_Parameters)
rownames(Data_Bio_Driven) <- NULL
write.csv(Data_Bio_Driven,"Biological Network with Data Driven Results fused.csv")
## Output : Data Driven Cluster Parameters added


## Input : Module 2 Clusters and TSEA Network
source("Module_4(Data Driven Network).R")
## Output : Network with Data Driven Clusters

