#Function to Incorporate clusters in the new data frame
#Input : data frame , list of clusters , print PCA Loadings used to combine the features
#Output : data frame as matrix
incorporate_groups <- function(df,group_list,printyn){
  if(printyn==1){
    print("Printing PCA Loadings used to combine the features")
  }
  grpnum = 1
  names_grp = c()
  pca_grp = c()
  for (i in group_list) {
    if(printyn==1){
      print(sprintf("Group %s : ",grpnum))
    }
    df_rest = df[, !(names(df) %in% i)]
    df_group = df[, (names(df) %in% i)]
    res.pca <- prcomp(df_group)
    pca.loadings <- res.pca$rotation
    col_name = c()
    new_col = 0
    k = 1
    for (j in colnames(df_group)) {
      new_col = new_col + (df[,j]*pca.loadings[k][1])
      names_grp = c(names_grp,j)
      pca_grp = c(pca_grp,pca.loadings[k][1])
      if(printyn==1){
        print(j)
        print(pca.loadings[k][1])}
      col_name = c(col_name,j)
      k = k + 1
    }
    col_name = paste(col_name, sep = "", collapse = "~")
    df_rest$col_name = new_col
    names(df_rest)[names(df_rest) == "col_name"] <- col_name
    df = df_rest
    grpnum = grpnum + 1 }
  grp_PCA <- cbind(names_grp,pca_grp)
  return( list(data.matrix(df),grp_PCA) )
}