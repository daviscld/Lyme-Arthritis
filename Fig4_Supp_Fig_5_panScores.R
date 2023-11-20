
#Z-score all together
X <- data_corrected
X <- scale(X, center = TRUE, scale = TRUE)

#Do NOT include ApoB100, which is an auto-antigen
IgG1_cols <- which(grepl("IgG1",colnames(X)) & !grepl("apoB", colnames(X)))
IgG2_cols <- which(grepl("IgG2",colnames(X)) & !grepl("apoB", colnames(X)))
IgG3_cols <- which(grepl("IgG3",colnames(X)) & !grepl("apoB", colnames(X)))
IgG4_cols <- which(grepl("IgG4",colnames(X)) & !grepl("apoB", colnames(X)))
IgA1_cols <- which(grepl("IgA1",colnames(X)) & !grepl("apoB", colnames(X)))
IgM_cols <- which(grepl("IgM",colnames(X)) & !grepl("apoB", colnames(X)))
secIgA_cols <- which(grepl("secIgA",colnames(X)) & !grepl("apoB", colnames(X)))
FcR2AH_cols <- which(grepl("FcR2AH",colnames(X)) & !grepl("apoB", colnames(X)))
FcR2B_cols <- which(grepl("FcR2B",colnames(X)) & !grepl("apoB", colnames(X)))
FcR3AF_cols <- which(grepl("FcR3AF",colnames(X)) & !grepl("apoB", colnames(X)))
FcR3B_cols <- which(grepl("FcR3B",colnames(X)) & !grepl("apoB", colnames(X)))
SNA_cols <- which(grepl("SNA",colnames(X)) & !grepl("apoB", colnames(X)))
RCA_cols <- which(grepl("RCA",colnames(X)) & !grepl("apoB", colnames(X)))
ADCD_cols <- which(grepl("ADCD",colnames(X)) & !grepl("apoB", colnames(X)))

##IGG PAN SCORES
#Create pan-antigen IgG1 score
panIgG1 <- as.data.frame(X[,IgG1_cols])
panIgG1$TotalG1 = rowSums(panIgG1)
#Create pan-antigen IgG2 score
panIgG2 <- as.data.frame(X[,IgG2_cols])
panIgG2$TotalG2 = rowSums(panIgG2)
#Create pan-antigen IgG3 score
panIgG3 <- as.data.frame(X[,IgG3_cols])
panIgG3$TotalG3 = rowSums(panIgG3)
#Create pan-antigen IgG4 score
panIgG4 <- as.data.frame(X[,IgG4_cols])
panIgG4$TotalG4 = rowSums(panIgG4)

#Create pan-antigen IgA1 score
panIgA1 <- as.data.frame(X[,IgA1_cols])
panIgA1$TotalA1 = rowSums(panIgA1)
#Create pan-antigen secIgA score
pansecIgA <- as.data.frame(X[,secIgA_cols])
pansecIgA$TotalsecIgA = rowSums(pansecIgA)

#Create pan-antigen IgM score
panIgM <- as.data.frame(X[,IgM_cols])
panIgM$TotalM = rowSums(panIgM)

#FC RECEPTOR SCORES
#Create pan-antigen FcR2AH score
panFcR2AH <- as.data.frame(X[,FcR2AH_cols])
panFcR2AH$TotalFcR2AH = rowSums(panFcR2AH)
#Create pan-antigen FcR2B score
panFcR2B <- as.data.frame(X[,FcR2B_cols])
panFcR2B$TotalFcR2B = rowSums(panFcR2B)
#Create pan-antigen FcR3AF score
panFcR3AF <- as.data.frame(X[,FcR3AF_cols])
panFcR3AF$TotalFcR3AF = rowSums(panFcR3AF)
#Create pan-antigen FcR3B score
panFcR3B <- as.data.frame(X[,FcR3B_cols])
panFcR3B$TotalFcR3B = rowSums(panFcR3B)

#FUNCTIONS SCORES
#Create pan-antigen SNA score
panSNA <- as.data.frame(X[,SNA_cols])
panSNA$TotalSNA = rowSums(panSNA)
#Create pan-antigen RCA score
panRCA <- as.data.frame(X[,RCA_cols])
panRCA$TotalRCA = rowSums(panRCA)
#Create pan-antigen ADCD score
panADCD <- as.data.frame(X[,ADCD_cols])
panADCD$TotalADCD = rowSums(panADCD)

panDF <- data.frame(panA1 = panIgA1$TotalA1,
                    panM = panIgM$TotalM,
                    panG1 = panIgG1$TotalG1,
                    panG2 = panIgG2$TotalG2,
                    panG3 = panIgG3$TotalG3,
                    #panG4 = panIgG4$TotalG4, #all zeroes
                    panFcR2AH = panFcR2AH$TotalFcR2AH,
                    panFcR2B = panFcR2B$TotalFcR2B,
                    panFcR3AF = panFcR3AF$TotalFcR3AF,
                    panFcR3B = panFcR3B$TotalFcR3B,
                    panSNA = panSNA$TotalSNA,
                    panRCA = panRCA$TotalRCA,
                    panADCD = panADCD$TotalADCD,
                    phenotype = phenotype,
                    fluid = fluid)

panDF_long <- reshape2::melt(cbind(as.data.frame(panDF),
                                   patientID = patientID_multilevel))

features <- unique(panDF_long$variable)

my_comparisons <- list(c("Joint Fluid", "Serum"))

for (phen in c("Responsive", "Refractory")){
  for (feature in features){

    if (grepl("ADCP", feature)){
      ylabel <- "sum of scaled MFI"
    } else {
      ylabel <- "sum of scaled log10 MFI"
    }

    testtest <- panDF_long[panDF_long$variable == feature & panDF_long$phenotype == phen,]
    plt <- ggplot(testtest, aes(x = fluid, y = as.numeric(value)), color = fluid) +
      geom_line(aes(group = patientID),
                color="grey",
                arrow = arrow(type = "closed",
                              length=unit(0.075, "inches"))) +
      geom_boxplot(aes(color = fluid, fill = fluid),
                   alpha = 0.5) +
      geom_point(aes(color = fluid), alpha = 1, size = 8) +
      scale_color_manual(values = c('#6666ae','#cb9ac6')) +
      scale_fill_manual(values = c('#6666ae','#cb9ac6')) +
      ylab(ylabel) +
      xlab(paste0(phen, " Fluid")) +
      labs(title = feature) +
      stat_compare_means(comparisons = my_comparisons,
                         paired = TRUE, #ONLY FOR JF v SERUM
                         symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                            symbols = c("****", "***", "**", "*", "ns")),
                         hide.ns = F) +
      theme_classic() +
      theme(legend.position = "none", axis.ticks.x=element_blank(),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 17),
            axis.text.x = element_text(size = 17),
            axis.title.x = element_text(size = 20),
            plot.title = element_text(hjust=0.5, size = 20))
    print(plt)
  }
}

#Supp Fig 5
panDF_long <- reshape2::melt(cbind(as.data.frame(panDF),
                                   patientID = patientID_multilevel))

my_comparisons <- list(c("Joint Fluid", "Serum"))
for (phen in c("Refractory", "Responsive")){

  testtest <- panDF_long[panDF_long$phenotype == phen,]
  plt <- ggplot(testtest, aes(x = fluid, y = as.numeric(value)), color = fluid) +
    geom_line(aes(group = patientID),
              color="grey",
              arrow = arrow(type = "closed",
                            length=unit(0.075, "inches"))) +
    geom_boxplot(aes(color = fluid, fill = fluid),
                 alpha = 0.5) +
    geom_point(aes(color = fluid), alpha = 1, size = 3) +
    scale_color_manual(values = c('#6666ae','#cb9ac6')) +
    scale_fill_manual(values = c('#6666ae','#cb9ac6')) +
    ylab("sum of scaled log10 MFI") +
    xlab(phen) +
    facet_wrap(.~variable, scales = "free_x", ncol = 4) +
    stat_compare_means(comparisons = my_comparisons,
                       paired = TRUE, #ONLY FOR JF v SERUM
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                          symbols = c("****", "***", "**", "*", "ns")),
                       hide.ns = F) +
    theme_classic() +
    theme(legend.position = "none", axis.ticks.x=element_blank(),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 17),
          axis.text.x = element_text(size = 8), #angle = 315, vjust = 0.1),
          plot.title = element_text(hjust=0.5, size = 20))
  print(plt)
}

#Corrected p values

my_comparisons <- list(c("Joint Fluid", "Serum"))
feature_list <- as.character(unique(panDF_long$variable))

#For responsive patients
p_paired_feats_resp <- list()
i = 0
for (feature in feature_list){
  i = i+1
  test <- panDF_long[which(grepl(feature, panDF_long$variable)),]
  test <- test[which(test$phenotype == "Responsive"),]
  test$group <- rep(NA, nrow(test))
  res <- wilcox.test(value ~ fluid, data = test, paired = TRUE)
  p_paired_feats_resp[i] <- res$p.value
}

p_paired_feats_resp <- unlist(p_paired_feats_resp)

pvalue_paired_BH_resp <- round(p.adjust(p_paired_feats_resp, "BH"),3)
indFeatures_paired <- which(pvalue_paired_BH_resp < 0.05 , pvalue_paired_BH_resp)
sig_paired_feat_resp <- feature_list[indFeatures_paired]
print("significant features: responsive jf vs serum")
print(sig_paired_feat_resp)
print("corrected p-values: will need to manually adjust on figures")
print(pvalue_paired_BH_resp[which(pvalue_paired_BH_resp < 0.05 , pvalue_paired_BH_resp)])

p_paired_feats_ref <- list()
i = 0

#For refractory patients
for (feature in feature_list){
  i = i+1
  test <- panDF_long[which(grepl(feature, panDF_long$variable)),]
  test <- test[which(test$phenotype == "Refractory"),]
  test$group <- rep(NA, nrow(test))
  res <- wilcox.test(value ~ fluid, data = test, paired = TRUE)
  p_paired_feats_ref[i] <- res$p.value
}

p_paired_feats_ref <- unlist(p_paired_feats_ref)

pvalue_paired_BH_ref <- round(p.adjust(p_paired_feats_ref, "BH"),3)
indFeatures_paired <- which(pvalue_paired_BH_ref < 0.05 , pvalue_paired_BH_ref)
sig_paired_feat_ref <- feature_list[indFeatures_paired]
print("significant features: refractory jf vs serum")
print(sig_paired_feat_ref)
print("corrected p-values: will need to manually adjust on figures")
print(pvalue_paired_BH_ref[which(pvalue_paired_BH_ref < 0.05 , pvalue_paired_BH_ref)])
