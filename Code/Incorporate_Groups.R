#Function to Incorporate clusters in the new data frame
#Input : data frame , list of clusters
#Output : data frame
incorporate_groups <- function(df,group_list){
  for (i in group_list) {
    if(length(i) > 1){
      df_rest = df[, !(names(df) %in% i)]
      df_group = df[, (names(df) %in% i)]
      res.pca <- prcomp(df_group)
      var_explain = summary(res.pca)[["importance"]][2][1]
      if(var_explain <0.80){print("Variance Explained by PC1 for the following group is < 80%")
        print(i) 
        print(var_explain)}
      pca.loadings <- res.pca$rotation
      col_name = c()
      new_col = 0
      k = 1
      for (j in colnames(df_group)) {
        new_col = new_col + (df[,j]*pca.loadings[k][1])
        #col_name = c(col_name,substr(j,1,5))
        col_name = c(col_name,j)
        k = k + 1
      }
      col_name = paste(col_name, sep = "", collapse = "~")
      df_rest$col_name = new_col
      names(df_rest)[names(df_rest) == "col_name"] <- col_name
      df = df_rest}
    
  }
  return(data.matrix(df))
}