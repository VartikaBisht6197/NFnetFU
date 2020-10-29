## Input : List of Microbes

print("----Microbiome Analyst----")

#MicrobiomeAnalyst
mbSet<-Init.mbSetObj()
mbSet<-SetModuleType(mbSet, "tsea")
mbSet<-Setup.MapData(mbSet, Microbes)
mbSet<-CrossReferencing(mbSet, "mixed")
mbSet<-SetTaxonSetLib(mbSet, "host_int")
mbSet<-CalculateHyperScore(mbSet)
network_disease_res <- current.msetlib
write.csv(network_disease_res,"Mix Taxa TSEA Results.csv")

print("Mix Taxa TSEA Results Calculated")

#Get the information for disease mentioned
index_OTU_res <- c()
for(i in 1:length(current.msetlib$name)){
  dis_check <- c()
  for (j in disease) {
    dis_check <- c(dis_check,grepl(j,current.msetlib$name[i], fixed = TRUE))
  }
  if(any(dis_check == TRUE)){
    index_OTU_res <- c(index_OTU_res,i)
  }
}

print("Mix Taxa TSEA Disease Specific Results Calculated")

#TSEA Results for particular disease
network_disease_res <- current.msetlib[index_OTU_res,]
write.csv(network_disease_res,"Mix Taxa TSEA Disease Specific Results.csv")


#Adjacency Matrix for Network

print("Calculating Adjacency Matrix for Network")

#Initialisation
adjacency_matrix <- matrix(data=0,nrow = length(Microbes), ncol = length(Microbes))
colnames(adjacency_matrix) <- Microbes
rownames(adjacency_matrix) <- Microbes

#Filling in the Values
for(i in 1:length(current.mset[index_OTU_res])){
  inter_OTU <- c()
  for(k in current.mset[index_OTU_res][i]){
    for (l in Microbes) {
      if(any(grepl(l,k, fixed = TRUE)) == TRUE){
        inter_OTU <- c(inter_OTU,l)
      }
    }
  }
  if(length(inter_OTU)>0){
    combini <- expand.grid(inter_OTU,inter_OTU)
    for (j in 1:dim(combini)[1]) {
      adjacency_matrix[as.character(combini[j,1]),as.character(combini[j,2])] <- adjacency_matrix[as.character(combini[j,1]),as.character(combini[j,2])] + 1
    }
  }
}

#Node size indicating frequency of occurance in the study
lit_wt <- diag(adjacency_matrix)
write.csv(lit_wt,"Node Size (TSEA).csv")

#Upper Triangular Adjacency Matrix
adjacency_matrix[upper.tri(adjacency_matrix, diag = TRUE)] <- 0
write.csv(adjacency_matrix,"Upper Triangular Adjacency Matrix (TSEA).csv")

#Making Edge Matrix for the Network
OTU_melt_matrix <- as.data.frame(melt(adjacency_matrix))
del <- c()
for (i in 1:length(OTU_melt_matrix$value)) {
  if(OTU_melt_matrix$value[i]==0)
    del <- c(del,i)
}
edge_matrix <- as.matrix(OTU_melt_matrix[-del,][,c(colnames(OTU_melt_matrix)[1],colnames(OTU_melt_matrix)[2])])
rownames(edge_matrix) <- seq(1,dim(edge_matrix)[1])

#Making the Network ( LAYOUT : layout_in_circle, layout_with_fr, layout_nicely(g) )
g <- graph_from_edgelist(edge_matrix, directed = FALSE)
E(g)$weight <- OTU_melt_matrix[-del,]$value
comps <- components(g)$membership
colbar <- colorRampPalette(brewer.pal(8, "Set3"))(max(comps)+1)
V(g)$color <- colbar[comps+1]
Legends <-cbind(1:length(vertex_attr(g)$name),vertex_attr(g)$name)
vertex_attr(g)$name <- 1:length(vertex_attr(g)$name)
vertex_attr(g)
vertex_wt <- c()
for (i in Legends[,2]) {
  vertex_wt <- c(vertex_wt,lit_wt[i])
}

#Legend with important information about the network
Network_Info <- cbind(Legends,as.numeric(vertex_wt))
colnames(Network_Info) <- c("Node","Microbe Names","Node Size")
write.csv(Network_Info,"Network Legends.csv")

vertex_wt <- vertex_wt + rep(10,length(vertex_wt))
tiff("Biological Network.tiff", width = 10, height = 10, units = 'in', res = 300)
plot(g, layout=layout_in_circle, vertex.size=vertex_wt,edge.width = E(g)$weight)
dev.off()


## Output : Network and Network Legends with Node size (Network_Info)
