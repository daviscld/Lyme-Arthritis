
#remember to run the correct data loading .R (Fig_2_5_Data_Load.R) and NOT the other loading .R prior to this!!!

#Run by antigen or by feature, with 4 flower plots (serum ref, serum res, JF ref, JF res)
#Z-score all together
X <- data_corrected
X <- scale(X, center = TRUE, scale = TRUE)
sep <- phenotype
sep <- droplevels(sep)
sep2 <- fluid
sep2 <- droplevels(sep2)

#RUN FOR BOTH BY ANTIGEN AND BY FEATURE
mypal <- colorRampPalette( brewer.pal( 9 , "Spectral" ) )

X_percent <- X
for (ind_feat in 1:ncol(X_percent)) {
  X_percent[, ind_feat] <- percent_rank(X_percent[, ind_feat])
}

X_percent_ref_JF <- X_percent[which(fluid == "Joint Fluid" & phenotype == "Refractory"),]
X_percent_ref_s <- X_percent[which(fluid == "Serum" & phenotype == "Refractory"),]
X_percent_res_JF <- X_percent[which(fluid == "Joint Fluid" & phenotype == "Responsive"),]
X_percent_res_s <- X_percent[which(fluid == "Serum" & phenotype == "Responsive"),]

X_percent_ref_JF <- melt(X_percent_ref_JF)
colnames(X_percent_ref_JF)<-c('sample','variable','value')
X_percent_ref_JF$variable <- factor(X_percent_ref_JF$variable, levels = unique(X_percent_ref_JF$variable))
X_percent_ref_JF$value <- as.numeric(X_percent_ref_JF$value)
X_percent_ref_JF<- aggregate(X_percent_ref_JF$value, by=list(X_percent_ref_JF$variable), FUN=mean)
colnames(X_percent_ref_JF)<-c('variable','value')
X_percent_ref_JF$variable <- gsub(" 100","", X_percent_ref_JF$variable)
X_percent_ref_JF$antigen <- unlist(strsplit(as.character(X_percent_ref_JF$variable), split=' '))[seq(0,2*nrow(X_percent_ref_JF),2)]
X_percent_ref_JF$value[is.na(X_percent_ref_JF$value)] <- 0

X_percent_ref_s <- melt(X_percent_ref_s)
colnames(X_percent_ref_s)<-c('sample','variable','value')
X_percent_ref_s$variable <- factor(X_percent_ref_s$variable, levels = unique(X_percent_ref_s$variable))
X_percent_ref_s$value <- as.numeric(X_percent_ref_s$value)
X_percent_ref_s<- aggregate(X_percent_ref_s$value, by=list(X_percent_ref_s$variable), FUN=mean)
colnames(X_percent_ref_s)<-c('variable','value')
X_percent_ref_s$variable <- gsub(" 100","", X_percent_ref_s$variable)
X_percent_ref_s$antigen <- unlist(strsplit(as.character(X_percent_ref_s$variable), split=' '))[seq(0,2*nrow(X_percent_ref_s),2)]
X_percent_ref_s$value[is.na(X_percent_ref_s$value)] <- 0

X_percent_res_JF <- melt(X_percent_res_JF)
colnames(X_percent_res_JF)<-c('sample','variable','value')
X_percent_res_JF$variable <- factor(X_percent_res_JF$variable, levels = unique(X_percent_res_JF$variable))
X_percent_res_JF$value <- as.numeric(X_percent_res_JF$value)
X_percent_res_JF<- aggregate(X_percent_res_JF$value, by=list(X_percent_res_JF$variable), FUN=mean)
colnames(X_percent_res_JF)<-c('variable','value')
X_percent_res_JF$variable <- gsub(" 100","", X_percent_res_JF$variable)
X_percent_res_JF$antigen <- unlist(strsplit(as.character(X_percent_res_JF$variable), split=' '))[seq(0,2*nrow(X_percent_res_JF),2)]
X_percent_res_JF$value[is.na(X_percent_res_JF$value)] <- 0

X_percent_res_s <- melt(X_percent_res_s)
colnames(X_percent_res_s)<-c('sample','variable','value')
X_percent_res_s$variable <- factor(X_percent_res_s$variable, levels = unique(X_percent_res_s$variable))
X_percent_res_s$value <- as.numeric(X_percent_res_s$value)
X_percent_res_s<- aggregate(X_percent_res_s$value, by=list(X_percent_res_s$variable), FUN=mean)
colnames(X_percent_res_s)<-c('variable','value')
X_percent_res_s$variable <- gsub(" 100","", X_percent_res_s$variable)
X_percent_res_s$antigen <- unlist(strsplit(as.character(X_percent_res_s$variable), split=' '))[seq(0,2*nrow(X_percent_res_s),2)]
X_percent_res_s$value[is.na(X_percent_res_s$value)] <- 0


##FLOWERS BY FEATURE
flower_feature_list <- c("ADCD","IgA1","secIgA","RCA","SNA", "IgG1", "IgG2",
                         "IgG3", "IgG4", "IgM", "FcR2AH", "FcR2B", "FcR3AF",
                         "FcR3B")
for (feat in flower_feature_list){
  X_percent_ref_JF_ft <- X_percent_ref_JF[which(grepl(feat,X_percent_ref_JF$variable)),]
  plt_JF_ref <- ggplot(X_percent_ref_JF_ft, aes(factor(antigen), value, fill = factor(antigen))) +
    coord_polar() +
    geom_hline(yintercept = seq(0,0.8, by = 0.2), color = 'grey', size = 1) +
    xlab(paste0(feat, "_", "JF_ref")) +
    geom_bar(stat='identity',width = 1)+
    scale_fill_manual( values = mypal(19) ) +
    #geom_vline(xintercept = Ig_breaks, color = 'grey', size = 1) +
    theme(panel.background = element_blank(),
          legend.position = "none")

  X_percent_ref_s_ft <- X_percent_ref_s[which(grepl(feat,X_percent_ref_s$variable)),]

  plt_serum_ref <- ggplot(X_percent_ref_s_ft, aes(factor(antigen), value, fill = factor(antigen))) +
    coord_polar() +
    geom_hline(yintercept = seq(0,0.8, by = 0.2), color = 'grey', size = 1) +
    xlab(paste0(feat, "_", "serum_ref")) +
    geom_bar(stat='identity',width = 1)+
    scale_fill_manual( values = mypal(19) ) +
    #geom_vline(xintercept = Ig_breaks, color = 'grey', size = 1) +
    theme(panel.background = element_blank(),
          legend.position = "none")

  X_percent_res_JF_ft <- X_percent_res_JF[which(grepl(feat,X_percent_res_JF$variable)),]

  plt_JF_res <- ggplot(X_percent_res_JF_ft, aes(factor(antigen), value, fill = factor(antigen))) +
    coord_polar() +
    geom_hline(yintercept = seq(0,0.8, by = 0.2), color = 'grey', size = 1) +
    xlab(paste0(feat, "_", "JF_res")) +
    geom_bar(stat='identity',width = 1)+
    scale_fill_manual( values = mypal(19) ) +
    #geom_vline(xintercept = Ig_breaks, color = 'grey', size = 1) +
    theme(panel.background = element_blank(),
          legend.position = "none")

  X_percent_res_s_ft <- X_percent_res_s[which(grepl(feat,X_percent_res_s$variable)),]

  plt_serum_res <- ggplot(X_percent_res_s_ft, aes(factor(antigen), value, fill = factor(antigen))) +
    coord_polar() +
    geom_hline(yintercept = seq(0,0.8, by = 0.2), color = 'grey', size = 1) +
    xlab(paste0(feat, "_", "serum_res")) +
    geom_bar(stat='identity',width = 1)+
    scale_fill_manual( values = mypal(19) ) +
    #geom_vline(xintercept = Ig_breaks, color = 'grey', size = 1) +
    theme(panel.background = element_blank(),
          legend.position = "none")

  grid.arrange(plt_JF_ref, plt_serum_ref, plt_JF_res, plt_serum_res)
}
