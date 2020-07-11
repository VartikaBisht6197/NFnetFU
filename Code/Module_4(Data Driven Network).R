## Input : Information about the clusters from Module 2 and TSEA Network

print("Calculating Edges and Nodes to be added for the Data Driven Network")

node_info <- c()
for (i in 1:length(groups_we_need)) {
  data_node <- c()
  for (j in 1:length(Cluster_OTU_name)) {
    counter <- any( groups_we_need[[i]] %in% Cluster_OTU_name[[j]] == TRUE )
    if(counter){
      data_node <- c(data_node,j)
    }
  }
  node_info <- append(node_info,list(data_node))
}


edge_no <- c() #Change Colour of Existing Edges
nodes_added <- c() #Add New Nodes To Network Info
node_no_added <- c() #Node Number to Network
edges_added <- c() #New Edges added
k = 1

for (i in node_info) {
  if( !is.null(i) ){
    #Existing Edges and Nodes
    if(length(i) > 1){
      for (j in 1:dim(edge_matrix)[1]) {
        if( any( edge_matrix[j,1] %in% Microbes[i] == TRUE ) ){
          edge_no <- c(edge_no,j)}}
    }else{ #New Edges and Nodes
      index <- which(as.character(Name_Change[groups_we_need[[k]],]) %in% Network_Info[,"Microbe Names"][i])
      new_node <- as.character(Name_Change[groups_we_need[[k]][-index],])
      if(length(new_node) > 0){
        nodes_added <- c(nodes_added,unique(new_node))
        for (l in unique(new_node)) {
          ap <- list(groups_we_need[[k]][which(as.character(Name_Change[groups_we_need[[k]],])==l)])
          Cluster_OTU_name <- append(Cluster_OTU_name,ap)
        }
        adm_nodes <- c(seq(dim(Network_Info)[1]+1,dim(Network_Info)[1]+length(new_node)),i)
        node_no_added <- c(node_no_added,adm_nodes[-length(adm_nodes)])
        adm <- combn(adm_nodes,2)
        
        for (m in 1:dim(adm)[2]) {
          edges_added <- c(edges_added,adm[1,m],adm[2,m])
        }
      }
    }
  }
  k = k + 1
}

#Making new Network using the original TSEA Network
gh <- g
col_edges <- rep("green",length(E(gh)))
col_edges[edge_no] <- "red"
gh <- add_vertices(gh,length(node_no_added),name=node_no_added)
gh <- add_edges(gh,edges_added)
col_edges <- c(col_edges,rep("red",length(edges_added)/2))
E(gh)$weight <- c(E(g)$weight,rep(1,(length(edges_added))/2))
vertex_wt_gh <- c(vertex_wt,rep(10,length(node_no_added)))

#Save the New Fused Network
tiff("Novel Data Driven Biological Network.tiff", width = 10, height = 10, units = 'in', res = 300)
plot(gh, layout=layout_in_circle, vertex.size=vertex_wt_gh,edge.width = E(gh)$weight,edge.color=col_edges)
dev.off()

print("Final Fused Network Saved!")

nodes_in_cluster <- c()
for(i in Cluster_OTU_name){
  nodes_in_cluster <- c(nodes_in_cluster,paste(i, collapse=', ' ))
}

Data_Bio_Driven_with_clusters <- as.matrix(Data_Bio_Driven)
if(!is.null(node_no_added)){
  for (i in 1:length(node_no_added)) {
    CP <- abs(feature_parameters[Cluster_OTU_name[node_no_added[i]][[1]],])
    entry_DB <- c(node_no_added[i],nodes_added[i],"1",CP)
    Data_Bio_Driven_with_clusters <- rbind(Data_Bio_Driven_with_clusters,entry_DB)
    rownames(Data_Bio_Driven_with_clusters) <- NULL
    temp_col <- c(colnames(Data_Bio_Driven_with_clusters) ,"Features in the Cluster")
    Data_Bio_Driven_with_clusters <- cbind(Data_Bio_Driven_with_clusters,nodes_in_cluster)
    colnames(Data_Bio_Driven_with_clusters) <- temp_col
    write.csv(Data_Bio_Driven_with_clusters,"Final Fused Network.csv")
  }
}


print("Final Fused Network Cluster Information Saved!")

## Output : Network with Data Driven Clusters (Data_Bio_Driven_with_clusters)
