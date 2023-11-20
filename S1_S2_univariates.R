#get longform data
patientID_testing <- gsub("jf|s", "", patientID)
dfBox <- reshape2::melt(cbind(as.data.frame(data_corrected),
                              phenotype = phenotype,
                              fluid = fluid,
                              patientID = patientID_testing))

#Calculate unmodified p-values using Wilcoxon tests (paired OR unpaired, depending)
feature_list <- as.character(unique(dfBox$variable))
np <- length(feature_list)

#FOR UNPAIRED COMPARISONS
#Ref vs res serum
pvalue_ref_v_res_serum <- rep(1, length=np)
subset_1 <-dfBox[dfBox$fluid == "Serum" & dfBox$phenotype == "Responsive",]
subset_2 <-dfBox[dfBox$fluid == "Serum" & dfBox$phenotype == "Refractory",]
subset <- rbind(subset_1, subset_2)

for (i in 1:np){
  featuredfeature <- feature_list[i]
  subset_sig <- subset[subset$variable == featuredfeature,]
  dfB <- wilcox.test(value ~ phenotype, data = subset_sig, paired = FALSE)
  pvalue_ref_v_res_serum[i] <- dfB[3]
}

#Ref vs res jf
pvalue_ref_v_res_jf <- rep(1, length=np)
subset_1 <-dfBox[dfBox$fluid == "Joint Fluid" & dfBox$phenotype == "Responsive",]
subset_2 <-dfBox[dfBox$fluid == "Joint Fluid" & dfBox$phenotype == "Refractory",]
subset <- rbind(subset_1, subset_2)

for (i in 1:np){
  featuredfeature <- feature_list[i]
  subset_sig <- subset[subset$variable == featuredfeature,]
  dfB <- wilcox.test(value ~ phenotype, data = subset_sig, paired = FALSE)
  pvalue_ref_v_res_jf[i] <- dfB[3]
}

#FOR PAIRED COMPARISONS
#ONLY take matched fluid for each phenotype
#Refractory jf vs s
pvalue_ref_jf_v_s <- rep(1, length=np)
subset_1 <-dfBox[dfBox$phenotype == "Refractory" & dfBox$fluid == "Joint Fluid",]
subset_2 <-dfBox[dfBox$phenotype == "Refractory" & dfBox$fluid == "Serum",]
subset_1 <- subset_1[which(subset_1$patientID %in% subset_2$patientID),]
subset_2 <- subset_2[which(subset_2$patientID %in% subset_1$patientID),]
subset <- rbind(subset_1, subset_2)

for (i in 1:np)
{
  featuredfeature <- feature_list[i]
  subset_sig <- subset[subset$variable == featuredfeature,]
  dfB <- wilcox.test(value ~ fluid, data = subset_sig, paired = TRUE)
  pvalue_ref_jf_v_s[i] <- dfB[3]
}
#Responsive jf vs s
pvalue_res_jf_v_s <- rep(1, length=np)
subset_1 <-dfBox[dfBox$phenotype == "Responsive" & dfBox$fluid == "Joint Fluid",]
subset_2 <-dfBox[dfBox$phenotype == "Responsive" & dfBox$fluid == "Serum",]
subset_1 <- subset_1[which(subset_1$patientID %in% subset_2$patientID),]
subset_2 <- subset_2[which(subset_2$patientID %in% subset_1$patientID),]
subset <- rbind(subset_1, subset_2)

for (i in 1:np)
{
  featuredfeature <- feature_list[i]
  subset_sig <- subset[subset$variable == featuredfeature,]
  dfB <- wilcox.test(value ~ fluid, data = subset_sig, paired = TRUE)
  pvalue_res_jf_v_s[i] <- dfB[3]
}

#Adjust p-values by Benjamini-Hochberg correction for multiple testing
pvalue_ref_v_res_serum_adj <- round(p.adjust(pvalue_ref_v_res_serum, "BH"),3)
pvalue_ref_v_res_jf_adj <- round(p.adjust(pvalue_ref_v_res_jf, "BH"),3)
pvalue_ref_jf_v_s_adj <- round(p.adjust(pvalue_ref_jf_v_s, "BH"),3)
pvalue_res_jf_v_s_adj <- round(p.adjust(pvalue_res_jf_v_s, "BH"),3)

#Which features survive significance testing?
indFeatures <- which(pvalue_ref_v_res_serum_adj < 0.05 , pvalue_ref_v_res_serum_adj)
sigfeat_ref_v_res_serum <- feature_list[indFeatures]
print("signficiant features refractory vs responsive serum:")
print(sigfeat_ref_v_res_serum)
print("p-values (corrected): manually change p-values on plots")
print(pvalue_ref_v_res_serum_adj[which(pvalue_ref_v_res_serum_adj < 0.05 , pvalue_ref_v_res_serum_adj)])

indFeatures <- which(pvalue_ref_v_res_jf_adj < 0.05 , pvalue_ref_v_res_jf_adj)
sigfeat_ref_v_res_jf <- feature_list[indFeatures]
print("signficiant features refractory vs responsive joint fluid:")
print(sigfeat_ref_v_res_jf)
print("p-values (corrected): manually change p-values on plots")
print(pvalue_ref_v_res_jf_adj[which(pvalue_ref_v_res_jf_adj < 0.05 , pvalue_ref_v_res_jf_adj)])

indFeatures <- which(pvalue_ref_jf_v_s_adj < 0.05 , pvalue_ref_jf_v_s_adj)
sigfeat_ref_jf_v_s <- feature_list[indFeatures]
print("signficiant features joint fluid vs serum, refractory:")
print(sigfeat_ref_jf_v_s)
print("p-values (corrected): manually change p-values on plots")
print(pvalue_ref_jf_v_s_adj[which(pvalue_ref_jf_v_s_adj < 0.05 , pvalue_ref_jf_v_s_adj)])

indFeatures <- which(pvalue_res_jf_v_s_adj < 0.05 , pvalue_res_jf_v_s_adj)
sigfeat_res_jf_v_s <- feature_list[indFeatures]
print("signficiant features joint fluid vs serum, responsive:")
print(sigfeat_res_jf_v_s)
print("p-values (corrected): manually change p-values on plots")
print(pvalue_res_jf_v_s_adj[which(pvalue_res_jf_v_s_adj < 0.05 , pvalue_res_jf_v_s_adj)])

#save these adjusted p-values- will use them to manually adjust ggplot's p-values later

#get patients with matched serum and jf samples
matched_ind <- list()
i = 0
for (patient in gsub("jf|s", "", patientID)){
  i = i+1
  if (length(patientID[which(grepl(patient, patientID))]) > 1){
    matched_ind[i] <- patient
  } else {
    matched_ind[i] <- NA
  }
}
keep <- unlist(matched_ind)[which(!is.na(unlist(matched_ind)))]

dfTest <- dfBox[which(dfBox$patientID %in% keep),]
matchedID <- gsub("jf|s", "", patientID)
matchedID <- matchedID[which(matchedID %in% keep)]

fluid_list <- c("serum", "jf")
for (fl in fluid_list){
  if (fl == "jf"){
    #for joint fluid - Fig 1, Supp. Fig 1
    subset_1 <-dfBox[dfBox$fluid == "Joint Fluid" & dfBox$phenotype == "Responsive",]
    subset_2 <-dfBox[dfBox$fluid == "Joint Fluid" & dfBox$phenotype == "Refractory",]
    subset <- rbind(subset_1, subset_2)
  } else {
    #for serum - Fig 1, Supp. Fig 2
    subset_1 <-dfBox[dfBox$fluid == "Serum" & dfBox$phenotype == "Responsive",]
    subset_2 <-dfBox[dfBox$fluid == "Serum" & dfBox$phenotype == "Refractory",]
    subset <- rbind(subset_1, subset_2)
  }

  dfTest <- subset
  dfTest$phenotype <- droplevels(dfTest$phenotype)

  #for subsetting by antigen
  dfTest$variable <- gsub("ApoB 100", "ApoB", dfTest$variable)
  dfTest$variable <- gsub("apoB100", "ApoB", dfTest$variable)

  dfTest$antigen <- unlist(str_split(dfTest$variable," "))[seq(
    2,length(unlist(str_split(dfTest$variable, " "))),2)]
  dfTest$feature <- unlist(str_split(dfTest$variable," "))[seq(
    1,length(unlist(str_split(dfTest$variable, " "))),2)]

  my_comparisons <- list(c("Refractory","Responsive"))

  feature_groups <- c("IgG1", "IgG2","IgG3","IgA1","IgM",
                      "ADCD","FcR2AH", "FcR2B", "FcR3AF",
                      "FcR3B", "SNA", "RCA")
  antigen_list <- unique(dfTest$antigen)

  for (feature in feature_groups){
    test <- dfTest[which(grepl(feature, dfTest$feature)),]
    featured_feature <- feature
    test$group <- rep(NA, nrow(test))
    test$group[which(test$phenotype == "Refractory")] <- "Refractory"
    test$group[which(test$phenotype == "Responsive")] <- "Responsive"

    if (grepl("ADCP", feature)){
      ylabel <- "MFI"
    } else {
      ylabel <- "log10 MFI"
    }

    # Statistical test

    plt <- ggplot(test, aes(x = group, y = as.numeric(value)), color = group) +
      geom_boxplot(aes(color = group, fill = group),
                   alpha = 0.5) +
      geom_point(aes(color = group), alpha = 1, size = 3) +
      scale_color_manual(values = c('#363795','#41bb93')) +
      scale_fill_manual(values = c('#363795','#41bb93')) +
      ylab(ylabel) +
      xlab(fl) +
      labs(title = feature) +
      facet_wrap(.~antigen, scales = "free_x", ncol = 6) +
      scale_y_continuous(
        breaks = c(0, 2, 4, 6),
        expand = expand_scale(mult = c(0, 0.05))
      ) +
      stat_compare_means(comparisons = my_comparisons,
                         paired = FALSE,
                         symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                            symbols = c("****", "***", "**", "*", "ns")),
                         hide.ns = T) +
      theme_classic() +
      theme(legend.position = "none", axis.ticks.x=element_blank(),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 17),
            axis.text.x = element_blank(), #angle = 315, vjust = 0.1),
            plot.title = element_text(hjust=0.5, size = 20))

    print(plt)
  }
}

#manually adjust p values based on those calculated &
#multiple hypothesis corrected previously
print("manually adjust p-values based on those calculated & multiple hypothesis corrected previously")
