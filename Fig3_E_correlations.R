#Spearman correlations between all features
tidyCors <- X_down %>%
  correlate(method = "spearman") %>%
  stretch()

tidyCors <- tidyCors[-which(is.na(tidyCors$r)),]

#Test correlation significance and store p value in new p variable
tidyCors$p <- 0
for (ind in 1:dim(tidyCors)[1]) {
  tmp <- cor.test(X_down[,tidyCors$x[ind]], X_down[,tidyCors$y[ind]],
                  method = "spearman", exact = FALSE)
  tidyCors$p[ind] <- tmp$p.value
}
#Multiple hypothesis correction via Benjamini-Hochberg
tidyCors$p <- p.adjust(tidyCors$p, method = "BH", n = length(tidyCors$p))

#Which features do you want to plot correlations to?
plotFeatures <- sel_features#$sel_features #selected features

#Pull out only significant (p<0.05) and highly correlated
#(abs(spearman corr) > 0.7) features
graphCors <- tidyCors %>%
  filter(p < 0.05 & (abs(r) > 0.7| is.na(r)) &
           (is.element(x, plotFeatures) | is.element(y, plotFeatures))) %>%
  graph_from_data_frame(directed = FALSE)

#Visualization options for network graphs
layout <- create_layout(graphCors, layout = 'igraph', algorithm = 'nicely')
#in_circle, nicely, with_kk are all good options
nodeColours <- vector(mode = "list", length = length(V(graphCors)$name))
# selected features colored by what they're highest in

mark <- rep(NA, length = length(V(graphCors)$name))
for (ind_feat in 1:length(V(graphCors)$name)) {
  tmp_mean <- rep(NA, length = nlevels(sep))
  for (ind_class in 1:nlevels(sep)) {
    tmp_mean[ind_class] <- mean(X_down[which(sep == levels(sep)[ind_class]),
                                       which(colnames(X_down) == V(graphCors)$name[ind_feat])])
  }
  mark[ind_feat] <- levels(sep)[which.max(tmp_mean)]
}
#mark  <- factor(mark, levels = levels(sep))


ind = 0
for (feat in mark){
  ind = ind+1
  if (feat == "Refractory"){
    mark[ind] <- '#363795'
  } else {
    mark[ind] <- '#41bb93'
  }
}

#selected feature nodes colored by group enrichment
nodeColours <- mark

# other features colored white
nodeColours[!is.element(V(graphCors)$name, plotFeatures)] <- 'white'

pltGraph <- ggraph(layout) +
  geom_edge_link(aes(color = r), edge_width = .5) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(colors = colorRampPalette(brewer.pal(n = 11,
                                                                   name = "Reds"))(100))+ #choose color palette for visualization
  geom_node_point(shape = 21, color = "gray",
                  fill = nodeColours, size = 10, stroke = 0.5) +
  geom_node_text(aes(label = name, fontface = "bold"), size = 3,
                 point.padding = NA,  box.padding = 0, force = 0.1, repel = T) +
  theme(aspect.ratio = 1) +
  theme_graph(background = "white", base_family = 'Helvetica')

plt <- plot(pltGraph)

print(plt)
