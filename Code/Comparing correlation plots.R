library(Hmisc)
library(jmuOutlier)
library(stringr)

D1.0 <- read.xlsx("NFnetFU_Dataset1_wih_labels.xlsx", sheetIndex = 1)[3:49]
D2.0 <- read.csv("/Users/Vartika_Bisht/NFnetFU/Data/NFnetFU_Dataset2_OTU_abundance(top 100).csv",header = TRUE, row.names = 1)
D3.0 <- read.csv("/Users/Vartika_Bisht/NFnetFU/Data/NFnetFU_Dataset3_OTU_abundance(top 100).csv",header = TRUE, row.names = 1)

D1.1 <- read.csv("Results D1/Rule Based matrix with Colinearity Handled.csv",header = TRUE, row.names = 1)
D2.1 <- read.csv("Results D2/Rule Based matrix with Colinearity Handled.csv",header = TRUE, row.names = 1)
D3.1 <- read.csv("Results D3/Rule Based matrix with Colinearity Handled.csv",header = TRUE, row.names = 1)

D1 <- read.csv("Results D1/Ruled Based Matrix.csv",header = TRUE, row.names = 1)
D2 <- read.csv("Results D2/Ruled Based Matrix.csv",header = TRUE, row.names = 1)
D3 <- read.csv("Results D3/Ruled Based Matrix.csv",header = TRUE, row.names = 1)

corrD1.0 <- cor(D1.0)
write.csv(corrD1.0,"D1 Correlation matrix original Dataset.csv")
corrD2.0 <- cor(D2.0)
write.csv(corrD2.0,"D2 Correlation matrix original Dataset.csv")
corrD3.0 <- cor(D3.0)
write.csv(corrD3.0,"D3 Correlation matrix original Dataset.csv")


corrD1 <- cor(D1)
write.csv(corrD1,"D1 Correlation matrix.csv")
corrD2 <- cor(D2)
write.csv(corrD2,"D2 Correlation matrix.csv")
corrD3 <- cor(D3)
write.csv(corrD3,"D3 Correlation matrix.csv")

p.mat.D1 <- rcorr(as.matrix(D1), type ="spearman")[["P"]]
write.csv(p.mat.D1,"D1 Correlation matrix p values.csv")
p.mat.D2 <- rcorr(as.matrix(D2), type ="spearman")[["P"]]
write.csv(p.mat.D2,"D2 Correlation matrix p values.csv")
p.mat.D3 <- rcorr(as.matrix(D3), type ="spearman")[["P"]]
write.csv(p.mat.D3,"D3 Correlation matrix p values.csv")

corrD1.1 <- cor(D1.1)
write.csv(corrD1.1,"D1 Correlation matrix after module 2.csv")
corrD2.1 <- cor(D2.1)
write.csv(corrD2.1,"D2 Correlation matrix after module 2.csv")
corrD3.1 <- cor(D3.1)
write.csv(corrD3.1,"D3 Correlation matrix after module 2.csv")

p.mat.D1.1 <- rcorr(as.matrix(D1.1), type ="spearman")[["P"]]
write.csv(p.mat.D1.1,"D1 Correlation matrix p values after module 2.csv")
p.mat.D2.1 <- rcorr(as.matrix(D2.1), type ="spearman")[["P"]]
write.csv(p.mat.D2.1,"D2 Correlation matrix p values after module 2.csv")
p.mat.D3.1 <- rcorr(as.matrix(D3.1), type ="spearman")[["P"]]
write.csv(p.mat.D3.1,"D3 Correlation matrix p values after module 2.csv")


norm(corrD1, type ="F") - norm(corrD1.0, type ="F")
norm(corrD2, type ="F") - norm(corrD2.0, type ="F")
norm(corrD3, type ="F") - norm(corrD3.0, type ="F")



comp_cor <- function(mat1, mat2){
  names_ref <- colnames(mat1)
  mat <- matrix(c("Variable_1","Variable_2","Value_in_original","Value_after_M1","Absolute_diff"),nrow = 1,ncol = 5)
  for (i in 1:(length(names_ref)-1)) {
    for(j in (i+1):length(names_ref)){
      if( (names_ref[i] %in% colnames(mat2)) && (names_ref[j] %in% colnames(mat2)) ){
      mat_entry <- c(names_ref[i],names_ref[j],mat1[names_ref[i],names_ref[j]],mat2[names_ref[i],names_ref[j]], abs(mat2[names_ref[i],names_ref[j]] - mat1[names_ref[i],names_ref[j]]) )
      mat <- rbind(mat,mat_entry)} else{
      mat_entry <- c(names_ref[i],names_ref[j],mat1[names_ref[i],names_ref[j]],"NA","NA")
      mat <- rbind(mat,mat_entry)
      }
    }
  }
  mat <- as.data.frame(mat, col_vector = mat[1,])
  mat <- mat[-1,]
  return(mat)
}

remove_points_in_the_end <- function(mat){
  col_n <- colnames(mat)
  row_n <- rownames(mat)
  col_last_chr <- str_sub(col_n,-1,-1)
  n_col_last_chr <- c()
  row_last_chr <- str_sub(row_n,-1,-1)
  n_row_last_chr <- c()
  
  for(i in 1:length(col_last_chr)){
    if(col_last_chr[i] == "."){
      n_col_last_chr <- c( n_col_last_chr, substr(col_n[i],1,nchar(col_n[i])-1))
    }else{
      n_col_last_chr <- c( n_col_last_chr, col_n[i])
    }
  }
  
  for(i in 1:length(row_last_chr)){
    if(row_last_chr[i] == "."){
      n_row_last_chr <- c( n_row_last_chr, substr(row_n[i],1,nchar(row_n[i])-1))
    }else{
      n_row_last_chr <- c( n_row_last_chr, row_n[i])
    }
  }
  
  colnames(mat) <- n_col_last_chr
  rownames(mat) <- n_row_last_chr
  
  return(mat)
}

comp1 <- comp_cor(remove_points_in_the_end(corrD1.0),remove_points_in_the_end(corrD1))
comp2 <- comp_cor(remove_points_in_the_end(corrD2.0),remove_points_in_the_end(corrD2))
comp3 <- comp_cor(remove_points_in_the_end(corrD3.0),remove_points_in_the_end(corrD3))

write.csv(comp1,"Comparison between before and after ANFIS (correlation plot) D1.csv")
write.csv(comp2,"Comparison between before and after ANFIS (correlation plot) D2.csv")
write.csv(comp3,"Comparison between before and after ANFIS (correlation plot) D3.csv")

