#For ease of plotting, features are manually selected here. These selected features
#could also be retrieved from the "sel_features" variable after running the appropriate PLSDA/mPLSDA

#features from mPLSDA separating responsive JF v serum - manual list
#"FcR3B p27" is left out because it contributes essentially only to the second PC
#and not the first, which is the basis of separation. This feature will be
#explicitly mentioned in the figure caption. It has a significant amount of
#correlations, mostly with other FcRs and RCA
res_fts <- c("IgA1 OspA", "secIgA B67", "IgG2 Flagellin", "IgG3 p39",
             "ADCD Flagellin", "ADCD OspC")
res_ft_data <- data_corrected[which(phenotype == "Responsive"),]
#-which(colMeans(data_corrected) < .01)]
responsive_patients <- patientID_multilevel[which(phenotype == "Responsive")]
res_ft_data <- multilevel_denoising(as.matrix(res_ft_data), as.factor(responsive_patients))
sep <- fluid[which(phenotype == "Responsive")]

tidyCors_res <- res_ft_data %>%
  correlate(method = "spearman") %>%
  stretch()

tidyCors_res <- tidyCors_res[-which(is.na(tidyCors_res$r)),]
#Test correlation significance and store p value in new p variable
tidyCors_res$p <- 0
for (ind in 1:dim(tidyCors_res)[1]) {
  tmp <- cor.test(res_ft_data[,tidyCors_res$x[ind]],
                  res_ft_data[,tidyCors_res$y[ind]],
                  method = "spearman", exact = FALSE)
  tidyCors_res$p[ind] <- tmp$p.value
}
#Multiple hypothesis correction via Benjamini-Hochberg
tidyCors_res$p <- p.adjust(tidyCors_res$p, method = "BH",
                           n = length(tidyCors_res$p))

#Which features do you want to plot correlations to?
plotFeatures <- res_fts

#Pull out only significant (p<0.05) and highly correlated
#(abs(spearman corr) > 0.7) features
graphCors_res <- tidyCors_res %>%
  filter(p < 0.01 & (abs(r) > 0.85 | is.na(r)) &
           (is.element(x, plotFeatures) | is.element(y, plotFeatures))) %>%
  graph_from_data_frame(directed = FALSE)

#Visualization options for network graphs
layout <- create_layout(graphCors_res, layout = 'igraph', algorithm = 'nicely')
#in_circle, nicely, with_kk are all good options
nodeColours <- vector(mode = "list", length = length(V(graphCors_res)$name))
# selected features colored by what they're highest in

mark <- rep(NA, length = length(V(graphCors_res)$name))
for (ind_feat in 1:length(V(graphCors_res)$name)) {
  tmp_mean <- rep(NA, length = nlevels(sep))
  for (ind_class in 1:nlevels(sep)) {
    tmp_mean[ind_class] <- mean(res_ft_data[which(sep == levels(sep)[ind_class]),
                                            which(colnames(res_ft_data) == V(graphCors_res)$name[ind_feat])])
  }
  mark[ind_feat] <- levels(sep)[which.max(tmp_mean)]
}
#mark  <- factor(mark, levels = levels(sep))

ind = 0
for (feat in mark){
  ind = ind+1
  if (feat == "Joint Fluid"){
    mark[ind] <- '#6666ae'
  } else {
    mark[ind] <- '#cb9ac6'
  }
}

#selected feature nodes colored by group enrichment
nodeColours <- mark

# other features colored white
nodeColours[!is.element(V(graphCors_res)$name, plotFeatures)] <- 'white'

pltGraph <- ggraph(layout) +
  geom_edge_link(aes(color = r), edge_width = .5) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11,
                                                                       name = "RdBu")))(100))+ #choose color palette for visualization
  geom_node_point(shape = 21, color = "gray",
                  fill = nodeColours, size = 5, stroke = 0.5, alpha = 1) +
  geom_node_text(aes(label = name), size = 5,
                 point.padding = NA,  box.padding = 0, force = 0.1, repel = T) +
  theme(aspect.ratio = 1) +
  theme_graph(background = "white", base_family = 'Helvetica')

plt <- plot(pltGraph)

print(plt)

#features from mPLSDA separating refractory JF v serum - manual list
ref_ft <- c("IgA1 DbpA", "secIgA B67", "IgG2 DbpA", "IgG2 Flagellin", "IgG2 p39", "FcR3B OspA",
            "FcR3B OspE", "ADCD Flagellin", "ADCD VlsE", "ADCD p39")
ref_ft_data <- data_corrected[which(phenotype == "Refractory"),]
responsive_patients <- patientID_multilevel[which(phenotype == "Refractory")]
ref_ft_data <- multilevel_denoising(as.matrix(ref_ft_data), as.factor(responsive_patients))
sep <- fluid[which(phenotype == "Refractory")]

tidyCors_ref <- ref_ft_data %>%
  correlate(method = "spearman") %>%
  stretch()

#Test correlation significance and store p value in new p variable
tidyCors_ref$p <- 0
for (ind in 1:dim(tidyCors_ref)[1]) {
  tmp <- cor.test(ref_ft_data[,tidyCors_ref$x[ind]],
                  ref_ft_data[,tidyCors_ref$y[ind]],
                  method = "spearman", exact = FALSE)
  tidyCors_ref$p[ind] <- tmp$p.value
}
#Multiple hypothesis correction via Benjamini-Hochberg
tidyCors_ref$p <- p.adjust(tidyCors_ref$p, method = "BH",
                           n = length(tidyCors_ref$p))

#Which features do you want to plot correlations to?
plotFeatures <- ref_ft

#Pull out only significant (p<0.05) and highly correlated
#(abs(spearman corr) > 0.7) features
graphCors_ref <- tidyCors_ref %>%
  filter(p < 0.01 & (abs(r) > 0.85 | is.na(r)) &
           (is.element(x, plotFeatures) | is.element(y, plotFeatures))) %>%
  graph_from_data_frame(directed = FALSE)

#Visualization options for network graphs
layout <- create_layout(graphCors_ref, layout = 'igraph', algorithm = 'in_circle')
#in_circle, nicely, with_kk are all good options
nodeColours <- vector(mode = "list", length = length(V(graphCors_ref)$name))
# selected features colored by what they're highest in

mark <- rep(NA, length = length(V(graphCors_ref)$name))
for (ind_feat in 1:length(V(graphCors_ref)$name)) {
  tmp_mean <- rep(NA, length = nlevels(sep))
  for (ind_class in 1:nlevels(sep)) {
    tmp_mean[ind_class] <- mean(ref_ft_data[which(sep == levels(sep)[ind_class]),
                                            which(colnames(ref_ft_data) == V(graphCors_ref)$name[ind_feat])])
  }
  mark[ind_feat] <- levels(sep)[which.max(tmp_mean)]
}
#mark  <- factor(mark, levels = levels(sep))

ind = 0
for (feat in mark){
  ind = ind+1
  if (feat == "Joint Fluid"){
    mark[ind] <- '#6666ae'
  } else {
    mark[ind] <- '#cb9ac6'
  }
}

#selected feature nodes colored by group enrichment
nodeColours <- mark

# other features colored white
nodeColours[!is.element(V(graphCors_ref)$name, plotFeatures)] <- 'white'

pltGraph <- ggraph(layout) +
  geom_edge_link(aes(color = r), edge_width = .5) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(colors = colorRampPalette(brewer.pal(n = 11,
                                                                   name = "Reds"))(100))+ #choose color palette for visualization
  geom_node_point(shape = 21, color = "gray",
                  fill = nodeColours, size = 7, stroke = 0.5) +
  geom_node_text(aes(label = name), size = 5,
                 point.padding = NA,  box.padding = 0, force = 0.1, repel = T) +
  theme(aspect.ratio = 1) +
  theme_graph(background = "white", base_family = 'Helvetica')

plt <- plot(pltGraph)

print(plt)
