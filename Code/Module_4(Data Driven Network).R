## Input : Information about the clusters from Module 2 and TSEA Network

print("Calculating Edges and Nodes to be added for the Data Driven Network")


edge_no <- c() #Change Colour of Existing Edges
nodes_added <- c() #Add New Nodes To Network Info
node_no_added <- c() #Node Number to Network
edges_added <- c() #New Edges added

if(length(groups_we_need)>0){
node_info <- c()
for (i in 1:length(groups_we_need)) {
  data_node <- c()
  for (j in 1:length(Cluster_OTU_name)) {
    counter <- any( groups_we_need[[i]] %in% Cluster_OTU_name[[j]] == TRUE )
    if(counter){
      ind <- which(rownames(Name_Change) == Cluster_OTU_name[[j]][1]) #can take anyone becasue Name_Change will be same
      indx <- which(as.character(Name_Change[ind,]) == Microbes)
      data_node <- c(data_node,indx)
    }
  }
  node_info <- append(node_info,list(data_node))
}


for (i in node_info) {
  if( !is.null(i) ){
    #Existing Edges and Nodes
    if(length(i) > 1){
      for (j in 1:dim(edge_matrix)[1]) {
        microbe_EDGE <- t(combn(Microbes[i], 2))
        for (w in nrow(microbe_EDGE)) {
          if( all( c(edge_matrix[j,1],edge_matrix[j,2])  == c(microbe_EDGE[w,1],microbe_EDGE[w,2]))||all( c(edge_matrix[j,1],edge_matrix[j,2])  == c(microbe_EDGE[w,2],microbe_EDGE[w,1])) ){
            edge_no <- c(edge_no,j)}
        }}
    }}}



in_net_info <- Network_Info[,"Microbe Names"]
for (k in 1:length(groups_we_need)) {
for (i in 1:nrow(Network_Info)) {
     #New Edges and Nodes
      index <- which(as.character(Name_Change[groups_we_need[[k]],]) %in% Network_Info[,"Microbe Names"][i])
      new_node <- as.character(Name_Change[groups_we_need[[k]][-index],])
      nadd <- as.character(groups_we_need[[k]][-index])
      new_node <- new_node[!(new_node %in% in_net_info)]
      nadd <- nadd[which(!(new_node %in% in_net_info) == TRUE)]
      if(length(new_node) > 1){
        nodes_added <- c(nodes_added,nadd)
        for (l in unique(new_node)) {
          ap <- list(groups_we_need[[k]][which(as.character(Name_Change[groups_we_need[[k]],])==l)])
          Cluster_OTU_name <- append(Cluster_OTU_name,ap)
        }
        
        adm_nodes <- c(seq(length(in_net_info)+1,length(in_net_info)+length(unique(new_node))))
        node_no_added <- c(node_no_added,adm_nodes)
        adm <- combn(adm_nodes,2)
        for (m in 1:dim(adm)[2]) {
          edges_added <- c(edges_added,adm[1,m],adm[2,m])
        }
        extra_ed <- as.character(Name_Change[groups_we_need[[k]],])
        extra_ed <- unique(extra_ed[(extra_ed %in% in_net_info)])
        extra_ed <- which( in_net_info %in% extra_ed)
        for (f in extra_ed) {
          for (fg in adm_nodes) {
            edges_added <- c(edges_added,f,fg)
          }
        }
        in_net_info <- c(in_net_info,unique(new_node))
      }else{
    if(length(new_node)==1){
      node_no_added <- c(node_no_added,length(in_net_info)+1)
      extra_ed <- as.character(Name_Change[groups_we_need[[k]],])
      extra_ed <- unique(extra_ed[(extra_ed %in% in_net_info)])
      extra_ed <- which( in_net_info %in% extra_ed)
      for (f in extra_ed) {
          edges_added <- c(edges_added,f,length(in_net_info)+1)
        }
    }
        in_net_info <- c(in_net_info,new_node)
        nodes_added <- c(nodes_added,nadd)
  }}}

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

write.csv(as.character(Cluster_OTU_name),"OTU-Microbes in each node.csv")

Data_Bio_Driven_with_clusters <- as.matrix(Data_Bio_Driven)
if(!is.null(node_no_added)){
  for (i in 1:length(node_no_added)) {
    CP <- abs(feature_scores[nodes_added[i],])
    entry_DB <- c(node_no_added[i],nodes_added[i],"1",CP)
    Data_Bio_Driven_with_clusters <- rbind(Data_Bio_Driven_with_clusters,entry_DB)
    rownames(Data_Bio_Driven_with_clusters) <- NULL
    write.csv(Data_Bio_Driven_with_clusters,"Final Fused Network.csv")
  }
}
}else{ 
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
  
  Data_Bio_Driven_with_clusters <- as.matrix(Data_Bio_Driven)
write.csv(Data_Bio_Driven_with_clusters,"Final Fused Network.csv")}
print("Final Fused Network Cluster Information Saved!")

## Output : Network with Data Driven Clusters (Data_Bio_Driven_with_clusters)

