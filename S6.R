
figures <- c("3B", "5A", "5B")

for (fig in figures){
  #get data for each of the 3 plots
  if (fig == "3B"){
    #Z-score only joint fluid data
    X_down <- data_corrected[which(fluid == "Joint Fluid"),]
    X_down <- as.matrix(scale(X_down, center = TRUE, scale = TRUE))
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    sep <- phenotype[which(fluid == "Joint Fluid")]
    sep <- droplevels(sep)
    opts = list(n_folds = 5, pt_trials = 0, rf_trials = 0, paired = F)
    val_X <- X_down

  } else if (fig == "5A"){
    #Z-score for matched PLSDA
    X_down <- data_corrected[which(phenotype == "Refractory"),]
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    refrac_patients <- patientID[which(phenotype == "Refractory")]
    refrac_patients <- gsub("s|jf","", refrac_patients)
    X_down <- multilevel_denoising(as.matrix(X_down), as.factor(refrac_patients))
    sep <- fluid[which(phenotype == "Refractory")]
    patientID_multilevel_ <- patientID_multilevel[which(phenotype == "Refractory")]

  } else if (fig == "5B"){
    #Z-score for matched PLSDA
    X_down <- data_corrected[which(phenotype == "Responsive"),]
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    responsive_patients <- patientID[which(phenotype == "Responsive")]
    responsive_patients <- gsub("s|jf","", responsive_patients)
    X_down <- multilevel_denoising(as.matrix(X_down), as.factor(responsive_patients))
    sep <- fluid[which(phenotype == "Responsive")]
    patientID_multilevel_ <- patientID_multilevel[which(phenotype == "Responsive")]
  }

  opts = list(n_folds = 5, pt_trials = 0, rf_trials = 0)
  y_numeric <- as.numeric(sep)
  y_numeric[which(y_numeric == 2)] <- 0
  n_trials <- 100 #if you change this value, you must also adjust the number
  #of plots going into the final ggroc manually (eg delete plotN - plotN+X)

  #Set options for feature selection approach, using LASSO feature selection
  opts_sel <- list(n_trials = 100, threshold = 0.8, return_count = FALSE)

  select <- function(X, y) { return(select_repeat(X, y,
                                                  selector = select_lasso,
                                                  options = opts_sel)) }
  method = list(select = select,
                train = train_ropls,
                predict = predict_ropls,
                score = score_accuracy)

  y_pred <- validate_repeat_plusAUCROC(as.matrix(X_down), y_numeric,
                                       method, opts, n_trials = n_trials)

  #Get dataframe of values from y_pred for each trial
  tmp <- as.data.frame(do.call(cbind, y_pred))

  y_pred_df <- data.frame(matrix(ncol = n_trials, nrow = nrow(X_down)))
  for (col in 1:ncol(tmp)){
    y_pred_df[,col] <- tmp[,col]
  }


  #plot ROC curves and add AUC value for mean
  test_roc <- roc(y_numeric ~ y_pred_df[,1], plot = TRUE, print.auc = FALSE, legacy.axes = T)

  #testing ggroc
  for (rep in 1:(n_trials-1)){
    rep <- rep + 1
    assign(paste("plot", rep, sep = ""),
           roc(y_numeric ~ y_pred_df[,rep],
               plot = TRUE, print.auc = FALSE, legacy.axes = T))
  }

  #add average curve (might have to manually add to plots)
  mean_roc <- roc(y_numeric ~ rowMeans(y_pred_df), plot = TRUE,
                  print.auc = TRUE, legacy.axes = T, col = "black")
  roc_list <- paste("plot", seq(2,n_trials), sep = "")
  n = 0
  for (roc in roc_list){
    n = n+1
    roc_list[n] <- get(roc)
  }

  p <- ggroc(list(test_roc, plot2,plot3,plot4,plot5,plot6,
                  plot7,plot8,plot9,plot10,plot11,plot12,plot13,plot14,
                  plot15,plot16,plot17,plot18,plot19,plot20,plot21,
                  plot22,plot23,plot24,plot25,plot26,plot27,plot28,plot29,plot30,
                  plot31,plot32,plot33,plot34,plot35,plot36,plot37,plot38,plot39,
                  plot40,plot41,plot42,plot43,plot44,plot45,plot46,plot47,plot48,
                  plot49,plot50,plot51,plot52,plot53,plot54,plot55,plot56,plot57,
                  plot58,plot59,plot60,plot61,plot62,plot63,plot64,plot65,plot66,
                  plot67,plot68,plot69,plot70,plot71,plot72,plot73,plot74,plot75,
                  plot76,plot77,plot78,plot79,plot80,plot81,plot82,plot83,plot84,
                  plot85,plot86,plot87,plot88,plot89,plot90,plot91,plot92,plot93,
                  plot94,plot95,plot95,plot96,plot97,plot98,plot99, plot100, mean_roc),
             legacy.axes = TRUE) +
    scale_color_manual(values = c(rep("deepskyblue",101),"black")) +
    xlab("FPR") + ylab("TPR") +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                 color="darkgrey", linetype="dashed") +
    theme_classic() +
    theme(legend.position = "none")

  print(p)
}
