install.packages("tnet")
library(tnet)
library(igraph)
library(sna)
setwd("/Users/yx921121/Documents/GWU/2018Spring/Social network analysis/Project")


#read edges csv
edge<-read.csv("member-edges.csv",header=TRUE, stringsAsFactors = FALSE, na.strings = "NA")
#remove default id number
edge<-edge[,2:4]
#set id as character since they represent members
edge[,1]<-as.character(edge[,1])
edge[,2]<-as.character(edge[,2])
#read nodes csv
nodes<-read.csv("meta-members.csv",header=TRUE, stringsAsFactors = FALSE, na.strings = "NA")
#set.seed(1234)
#nodes<- nodes[sample(nrow(nodes),500),]
#only analyze vertexes of edges in nodes.csv
edgenew=edge[edge$member1 %in% nodes$member_id,]
edgenew2=edge[edge$member2 %in% nodes$member_id,]

#get undirected full model
full<-graph_from_data_frame(d=edgenew2, directed = TRUE)
V(full)
full.adj<-get.adjacency(full)
full.matrix<-as.matrix(full.adj)
summary(full)

#since original data is too huge to plot, we will only use first 1,000 rows data to plot
partialplot1<-graph_from_data_frame(d=edgenew2[1:1000,], directed = TRUE)
partialplot1sim <- simplify(partialplot1, remove.multiple = F, remove.loops = T) 
plot.igraph(partialplot1sim,vertex.label.cex=0.5,vertex.size=5,vertex.label=NA,edge.arrow.size=.4)

dg <- decompose.graph(partialplot1) # returns a list of three graphs
plot(dg[[1]]) # plot e.g. the 1st one


# Compute the indegree and outdegree centrality for each node
deg_in <- igraph::degree(full, mode="in") 
deg_out <- igraph::degree(full, mode="out") 
#degree_net<- igraph::degree(full)

# Compute shortest paths between each pair of nodes. 
sp_in <- igraph::shortest.paths(full, mode='in')
sp_out <- igraph::shortest.paths(full, mode='out')
# Closeness centrality

closeness_in <- igraph::closeness(full, mode='in')
closeness_out <- igraph::closeness(full, mode='out')
#closeness_net <- igraph::closeness(full)
# Betweenness centrality measures the number of shortest paths
# going through a specific vertex; it is returned by the 
# betweenness() function.
between_net <- igraph::betweenness(full)
#between_net

# Eigenvector centrality gives greater weight to a node the more 
# it is connected to other highly connected nodes. A node
# connected to five high-scoring nodes will have higher 
# eigenvector centrality than a node connected to five low-scoring
# nodes. Thus, it is often interpreted as measuring a node's
# network importance.
eigenv_net <- igraph::evcent(full)
#eigenv_net
eigenv_net_vector <- igraph::evcent(full,directed=FALSE)$vector

central_net <- data.frame(V(full)$name, deg_in,deg_out,closeness_in,closeness_out, between_net, eigenv_net_vector)
head(central_net)

central_net[order(-central_net$deg_in),] 
central_net[order(-central_net$deg_out),] 
central_net[order(-central_net$closeness_in),] 
central_net[order(-central_net$closeness_out),]   
central_net[order(-central_net$between_net),]
central_net[order(-central_net$eigenv_net_vector),]
barplot(central_net$deg_in,names.arg=central_net$V.full..name)
boxplot(central_net$deg_in,main="distribution of deg_in")
barplot(central_net$deg_out,names.arg=central_net$V.full..name)
barplot(central_net$closeness_in,names.arg=central_net$V.full..name)
barplot(central_net$closeness_out,names.arg=central_net$V.full..name)
barplot(central_net$between_net,names.arg=central_net$V.full..name)
barplot(central_net$eigenv_net_vector,names.arg=central_net$V.full..name)

central_net[with(central_net, order(-deg_in, -deg_out, closeness_in, closeness_out, -between_net, -eigenv_net_vector)),]


#CORRELATIONS BETWEEN CENTRALITY MEASURES
cor(central_net[,2:7])

install.packages("corrplot")
library(corrplot)
corrplot(cor(central_net[,2:7]), order = "hclust", 
         tl.col = "black", tl.srt = 45)

#merge data and then used for glm
require(plyr)
home <- nodes[,c(1,3)]
colnames(home) <-c('V.full..name','city')
mergedata <- merge(x=central_net, y=home)
between.glm <- glm(between_net~deg_in+deg_out+closeness_in+closeness_out+eigenv_net_vector+city, data=mergedata)
summary(between.glm)
cor(mergedata[,2:8])

## Reachability
#Reachability can only be computed on one vertex at a time. To get graph-wide statistics, 
#change the value of "vertex" manually or write a for loop. (Remember that, unlike R objects, 
#igraph objects are numbered from 0.)
reachability <- function(g, m) {   
  reach_mat = matrix(nrow = vcount(partialplot1),                       
                     ncol = vcount(partialplot1))   
  for(i in 1:vcount(partialplot1)) {     
    reach_mat[i,] = 0     
    this_node_reach <- subcomponent(g, i, mode = m) # used "i" instead of "(i - 1)"      
    for(j in 1:(length(this_node_reach))) {       
      alter = this_node_reach[j] # removed "+ 1"       
      reach_mat[i, alter] = 1     
    }   
  }   
  return(reach_mat) 
} 

reach_full_in <- reachability(partialplot1, 'in')
reach_full_out <- reachability(partialplot1, 'out')

#reach_full_in
#reach_full_out
mean(reach_full_in)
mean(reach_full_out)
sd(reach_full_in)
sd(reach_full_out)

## Distance
# Compute shortest paths between each pair of nodes. 
sp_full_in <- shortest.paths(partialplot1, mode='in')
sp_full_out <- shortest.paths(partialplot1, mode='out')
sp_full_in[is.infinite(sp_full_in)]<-0
sp_full_out[is.infinite(sp_full_out)]<-0


#sp_full_in
#sp_full_out
mean(sp_full_in)
mean(sp_full_out)
sd(sp_full_in)
sd(sp_full_out)

# Density 
graph.density(partialplot1)

# Reciprocity
reciprocity(partialplot1)

# Transitivity
transitivity(partialplot1)

## COMMUNITY DETECTION
#fullundirect<-as.undirected(full,mode = "collapse")
#edges<-edge.betweenness.community(fullundirect)
#edges
#plot(as.dendrogram(edges))


# Function to layout  coreness
CorenessLayout <- function(g) {
  coreness <- graph.coreness(g);
  xy <- array(NA, dim=c(length(coreness), 2));
  
  shells <- sort(unique(coreness));
  for(shell in shells) {
    v <- 1 - ((shell-1) / max(shells));
    nodes_in_shell <- sum(coreness==shell);
    angles <- seq(0,360,(360/nodes_in_shell));
    angles <- angles[-length(angles)]; # remove last element
    xy[coreness==shell, 1] <- sin(angles) * v;
    xy[coreness==shell, 2] <- cos(angles) * v;
  }
  return(xy);
}

# compute coreness
coreness <- graph.coreness(full);
# assign colors
colbar <- rainbow(max(coreness));
# create layout
ll <- CorenessLayout(full);
# plot
plot(full, layout=ll, 
     vertex.size=5, 
     vertex.color=colbar[coreness], 
     vertex.frame.color=colbar[coreness],
     vertex.label=NA,
     main='Coreness')
