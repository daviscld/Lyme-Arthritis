if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("USGS-R/inlmisc", dependencies = TRUE)
library(inlmisc)

#Z-score all together
X <- data_corrected
X <- scale(X, center = TRUE, scale = TRUE)
sep <- phenotype
sep <- droplevels(sep)
sep2 <- fluid
sep2 <- droplevels(sep2)

X_pheat <- cbind(as.data.frame(X), fluid = sep2, phenotype = sep)
X_pheat <- X_pheat[order(fluid),]

#Remove IgG4, C1q, ADCP, ADCD from heatmap (just do titers and lectins)
X_pheat <- X_pheat[,-which(grepl("IgG4|C1q|ADCP|ADCD", colnames(X_pheat)))]

fluid_pheat<- X_pheat$fluid
pheno_pheat <- X_pheat$phenotype
X_pheat <- X_pheat[,1:(ncol(X_pheat)-2)] #remove fluid/phenotype columns

annotation_row = data.frame(fluid = fluid_pheat,
                            phenotype = pheno_pheat)
rownames(annotation_row) <- rownames(X_pheat)

col_X <- colnames(X_pheat)

feature_names <- unlist(str_split(col_X," "))[seq(
  1,length(unlist(str_split(col_X, " "))),2)]
antigen_names <- unlist(str_split(col_X," "))[seq(
  2,length(unlist(str_split(col_X, " "))),2)]
annotation_col <- data.frame(feature = feature_names,
                             antigen = antigen_names)
rownames(annotation_col) <- colnames(X_pheat)


my_colors <- list(
  fluid = c("Joint Fluid" = '#6666ae', "Serum" = '#cb9ac6'),
  phenotype = c("Refractory" = '#363795', "Responsive" = '#41bb93'),
  feature = magma(length(unique(feature_names)))[sample(1:12,12)],
  antigen = inlmisc::GetTolColors(18, scheme = "discrete rainbow"))
names(my_colors$feature) <- unique(feature_names)
names(my_colors$antigen) <- unique(antigen_names)

X_pheat[is.na(X_pheat)] <- as.double("NA")

#Quantile breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs,
                     probs = seq(0, 1, length.out = n),
                     na.rm = T)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(-X_pheat, n = 11)

#Column breaks
col_gaps <- c(0)
for (feat in feature_names){
  break_ind <- match(feat,feature_names) - 1
  col_gaps <- c(col_gaps, break_ind)
}

#Create heatmap of the z-scored data indicating group
p_zscore <- pheatmap(X_pheat, annotation_colors = my_colors,
                     annotation_row = annotation_row,
                     annotation_col = annotation_col,
                     gaps_col = unique(col_gaps),
                     gaps_row = length(fluid_pheat[which(fluid_pheat == "Joint Fluid")]),
                     #color = viridis(length(mat_breaks) - 1),
                     color = RColorBrewer::brewer.pal((length(mat_breaks) - 1), "BrBG"),
                     cluster_cols = F, cluster_rows = F,
                     breaks = mat_breaks,
                     cellwidth = 3, cellheight = 3, fontsize = 6,
                     border_color = "NA", treeheight_row = 0,
                     treeheight_col = 0)

print(p_zscore)
