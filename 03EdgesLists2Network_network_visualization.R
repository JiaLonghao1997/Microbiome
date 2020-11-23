library("igraph")

setwd("F:\\Zhaolab2020\\BioNetwork\\MicrobiomeNetwork\\MNDnetwork\\PRJEB17784\\co-network")

dd <- read.csv("03merge_sparcc_spiec_edge_lists.csv", header=TRUE)
print(dd[1:5,])
dd$color=NA
for(i in 1:nrow(dd)){
  if(dd$coefficient[i]>0){
    dd$color[i]=1
  }
  else{
    dd$color[i]=-1
  }
}
print(dd[1:5, ])
gg <- graph.data.frame(dd, directed=FALSE)
E(gg)
V(gg)
?graph.data.frame
clr <- as.factor(dd$color)
print(clr)
levels(clr) <- c("green", "red")
?levels
print(clr)
vertex_attr(gg)
edge_attr(gg)
print(E(gg))
node_color <- clr[1:gorder(gg)]
print(node_color)
E(gg)$color <- as.character(clr)
print(edge_attr(gg))
#V(gg)$color <- as.character(node_color)
print(vertex_attr(gg))
print(V(gg)$color)
V(gg)$size = degree(gg)*0.5
#V(gg)$color <- "black"
# igraph static plot
# plot(g, layout = layout.circle, vertex.label=NA)
?edgebundle
edgebundle(gg, padding=250, fontsize=14,tension=0.6, width=1000, ver)
#circle_plot <- edgebundle(gg, padding=250, fontsize=14,tension=0.6, width=1000)
#saveEdgebundle(circle_plot, "03Circle_plot_Control.html", selfcontained = TRUE)

library(plotrix)

# first we create a function which defines what happens at each iteration of our algorithm
edge_betweenness_step = function(mat){
  # graph the matrix
  net <- graph.adjacency(mat)
  # evaluate edge betweenness
  edgebets <- edge_betweenness(net)
  # identify the edge that has the highest betweenness
  to_delete <- which.max(edgebets)
  # delete that edge
  net_new <- delete.edges(net, to_delete)
  # evaluate the number of components
  net_new <- decompose(net_new)
  # if the number of components grew to be greater than 1, return the matrices of those two components
  if(length(net_new) > 1){
    return(lapply(net_new, get.adjacency))
  } else {
    # otherwise return the matrix that we started with, minus the deleted edge
    return(get.adjacency(net_new[[1]]))
  }
}

# next we will repeatedly run this step algorithm a certain number of times, at each step applying it to every component we have created
# so at first, we start with the whole graph, we delete the highest between edges until the network is split into two components. 
# then we apply the above to each of the network subsets
# we do this until we have made N effective cuts
edge_betweenness_clustering <- function(net, desired_depth = 15){
  # required packages
  require(reshape2)
  require(plotrix)
  # turn our net into a matrix
  net_mat <- get.adjacency(net)
  # run the first iteration of edge_betweenness and save the result in a list
  net_temp <- list(edge_betweenness_step(net_mat))
  # while actual depth is < desired_depth
  # depth signifies how many effective cuts have been made
  # where an effective cut is any cut that divides a component into at least two new components
  while(listDepth(net_temp) < desired_depth){
    # apply edge_betweenness recursively to every component in the list
    net_temp <- rapply(net_temp, edge_betweenness_step, how = "list")
  }
  #get the row.names of every matrix (i.e. the actors who are in each component/group at the end of the iterations)
  groups <- lapply(unlist(net_temp), row.names)
  # name the groups according to their order
  names(groups) <- as.character(1:length(groups))
  # use melt to produce a person to group data.frame
  memberships <- melt(groups)
  # convert value (i.e. id) to a character vector from factor
  memberships$value <- as.character(memberships$value)
  # reorder the memberships data.frame so that it matches the order of vertices in the original network
  memberships <- memberships[match(V(net)$name, memberships$value),]
  # construct the communities object using this helpful function provided by igraph
  output <- make_clusters(net, 
                          membership = as.numeric(memberships$L1), 
                          algorithm = "edge_betweenness_attempt", 
                          modularity = TRUE)
  # return the communities object
  return(output)
}

plot(gg, 
     vertex.label = NA, 
     vertex.size=6, 
     edge.arrow.size = .1, 
     mark.groups = edge_betweenness_clustering(gg, 3))
# plot(gg)
# 
# 
# plot(gg, edge.arrow.size=.2, edge.curved=0,
#      
#      vertex.color="orange", vertex.frame.color="#555555",
#      
#      vertex.label.color="black",
#      
#      vertex.label.cex=.7) 
# V(gg)
# E(gg)
