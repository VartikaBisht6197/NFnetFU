if(data_analysis =="Dataset_1"){
  # Load Dataset 1
  Data_Set_1 <- read.csv("Data_Set_1.csv")
  # Choose Microbiome Data
  df_data1 <- Data_Set_1[3:22]
  data1 <- data.matrix(Data_Set_1[3:22])
  # Create labels for prediction ( 2nd column )
  label_dat <- label_creator(Data_Set_1)}


if(data_analysis =="Dataset_2"){
  # Load Dataset 2
  Data_Set_1 <- read.xlsx("Data_all_micobiome.xlsx",1)
  # Choose Microbiome Data
  df_data1 <- Data_Set_1[3:49]
  data1 <- data.matrix(Data_Set_1[3:49])
  # Create labels for prediction ( 2nd column )
  label_dat <- label_creator(Data_Set_1)}


if(data_analysis =="Dataset_3"){
  # Load Dataset 3
  Data_Set_1 <- read.table(file = "Data_set_git_3.tsv", sep = '\t', header = TRUE)
  #Choose Microbiome Data
  data1 <- data.matrix(Data_Set_1[2:5207])
  df_data1 <- Data_Set_1[2:5207]
  #Create labels for prediction ( 2nd column )
  label_dat <- Data_Set_1[1]}


if(data_analysis =="Stimulated_Dataset"){
  source("Simulated_data_set.R")
  # Stimulated Dataset
  stimulated_data <- data[1:20,2:81]
  #Choose Microbiome Data
  data1 <- data.matrix(stimulated_data)
  df_data1 <- stimulated_data
  #Create labels for prediction ( 2nd column )
  label_dat <- as.numeric(class_y[1:20])}

if(data_analysis =="Dataset_adenoma"){
  # Load Dataset 1
  Data_Set_1 <- read.table("final.0.03.subsample.csv", header = TRUE)
  meta_data <- read.csv("meta date.csv")
  index = c()
  for (i in 1:length(meta_data$sample)) {
    index_need <- match(meta_data$sample[i],Data_Set_1$Group)
    index <- c(index,index_need)
  }
  Data_Set_1 <- Data_Set_1[index,]
  # Choose Microbiome Data
  df_data1 <- Data_Set_1[4:6393]
  # Create labels for prediction ( 2nd column )
  label_dat <- as.numeric(factor(meta_data$dx))
  results <- rfe(df_data1, meta_data$dx, rfeControl=control)
  RFFeatures <- results[["optVariables"]][1:100]
  index = c()
  for (i in 1:length(RFFeatures)) {
    index_need <- match(RFFeatures[i],names(df_data1))
    index <- c(index,index_need)
  }
  df_data1 <- df_data1[,index]
  
  data1 <- data.matrix(df_data1)
}