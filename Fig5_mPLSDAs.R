#This will return Fig5A, B, C, D as well as cross-validation scores (may take a while to run)

opts_plot1 <- list(df_features = df_id,
                   loading_alpha = 1, # transparency for the loadings
                   score_alpha = 1, # transparency for the scores
                   color_features = "fluid",
                   colors = my_colors,
                   y_name = "fluid")

opts_plot3 <- list(df_features = df_features,
                   loading_alpha = 1, # transparency for the loadings
                   score_alpha = 1, # transparency for the scores
                   LV_ind = c(1,2), # which LVs to plot
                   n_LV = 2,
                   colors = my_colors,
                   y_name = "fluid")

#Remove features with too much missing data (>50% are 0) so they aren't chosen
#as features for PLSDA
remove_ind <- list()
n = 0
for (column in colnames(data_corrected)){
  n = n+1
  values <- data_corrected[,n]
  if (length(values[which(values < 0.2)]) > 42)
    remove_ind <- append(remove_ind, colnames(data_corrected)[n])
}
remove_ind <- unlist(remove_ind)

figure <- c("5A","5B")

for (fig in figure){
  if (fig == "5A"){
    #Z-score for matched PLSDA
    X_down <- data_corrected[which(phenotype == "Refractory"),]
    X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    refrac_patients <- patientID[which(phenotype == "Refractory")]
    refrac_patients <- gsub("s|jf","", refrac_patients)
    #multilevel denoising is what makes this a MATCHED PLS-DA, as the mean
    #of a single subject's two measurements is subtracted from both
    X_down <- multilevel_denoising(as.matrix(X_down), as.factor(refrac_patients))
    X_down <- as.data.frame(X_down)[-which(is.na(colMeans(X)))]
    sep <- fluid[which(phenotype == "Refractory")]
    patientID_multilevel_ref <- patientID_multilevel[which(phenotype == "Refractory")]
  } else if (fig == "5B"){
    #Z-score for matched PLSDA
    X_down <- data_corrected[which(phenotype == "Responsive"),]
    X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    responsive_patients <- patientID[which(phenotype == "Responsive")]
    responsive_patients <- gsub("s|jf","", responsive_patients)
    #multilevel denoising is what makes this a MATCHED PLS-DA, as the mean
    #of a single subject's two measurements is subtracted from both
    X_down <- multilevel_denoising(as.matrix(X_down), as.factor(responsive_patients))
    sep <- fluid[which(phenotype == "Responsive")]
    patientID_multilevel_res <- patientID_multilevel[which(phenotype == "Responsive")]
  }


  #Take features selected in 80% of 100 rounds of LASSO
  opts_sel <- list(n_trials = 100, threshold = 0.8, return_count = FALSE)


  sel_features <- select_repeat(as.matrix(X_down), sep,
                                selector = select_lasso,
                                options = opts_sel)

  X_sel <- X_down[, sel_features]


  # Perform an mPLS-DA using the selected features and plot the scores and loadings
  # Check number of latent variables and increase to 2 if <2 (for visualization purposes)

  opts_model <- list(n_LV = 2)
  model <- train_ropls(X_sel, sep, options = opts_model)
  ropls::getSummaryDF(model)
  plt_scores <- visualize_ropls_scores(model, sep, options = opts_plot1)
  print(plt_scores)

  # set additional options required to color code enrichment in the bar plot of the loadings
  opts_plot3$X <- X_sel
  opts_plot3$y <- sep
  opts_plot3$LV_ind <- 1
  opts_plot3$mark_enrichment <- TRUE
  plt_loadings_bar1 <- visualize_ropls_loadings_bar(model, options = opts_plot3)
  print(plt_loadings_bar1)
}

#cross validation
#Validation, Permutation Testing on Selected Features

for (fig in figure){
  if (fig == "5A"){
    #Z-score for matched PLSDA
    X_down <- data_corrected[which(phenotype == "Refractory"),]
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    refrac_patients <- patientID[which(phenotype == "Refractory")]
    refrac_patients <- gsub("s|jf","", refrac_patients)
    #multilevel denoising is what makes this a MATCHED PLS-DA, as the mean
    #of a single subject's two measurements is subtracted from both
    X_down <- multilevel_denoising(as.matrix(X_down), as.factor(refrac_patients))
    sep <- fluid[which(phenotype == "Refractory")]
    patientID_multilevel_ <- patientID_multilevel[which(phenotype == "Refractory")]
  } else if (fig == "5B"){
    #Z-score for matched PLSDA
    X_down <- data_corrected[which(phenotype == "Responsive"),]
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    responsive_patients <- patientID[which(phenotype == "Responsive")]
    responsive_patients <- gsub("s|jf","", responsive_patients)
    #multilevel denoising is what makes this a MATCHED PLS-DA, as the mean
    #of a single subject's two measurements is subtracted from both
    X_down <- multilevel_denoising(as.matrix(X_down), as.factor(responsive_patients))
    sep <- fluid[which(phenotype == "Responsive")]
    patientID_multilevel_ <- patientID_multilevel[which(phenotype == "Responsive")]
  }

  #To use paired cross-validation for mPLSDA, you MUST sort dataframe by sample type
  #as you need to have all of one sample type followed by all of the paired type
  val_X <- as.data.frame(X_down)
  val_X$fluid <- sep
  val_X$ID <- patientID_multilevel_
  val_X <- val_X[order(val_X$fluid),]
  val_y <- sort(sep)

  val_ID <- as.data.frame(patientID_multilevel_)
  val_ID$fluid <- sep
  val_ID <- val_ID[order(val_ID$fluid),]
  val_X <- val_X[,-c(ncol(val_X), (ncol(val_X)-1))]
  val_ID <- as.factor(val_ID[,-c(ncol(val_ID))])


  #Set options for feature selection approach, using LASSO feature selection
  opts_sel <- list(n_trials = 100, threshold = 0.8, return_count = FALSE)

  select <- function(X, y) { return(select_repeat(X, y,
                                                  selector = select_lasso,
                                                  options = opts_sel)) }
  method = list(select = select,
                train = train_ropls,
                predict = predict_ropls,
                score = score_accuracy)

  #How many total trials, permutation (pt) trials, and random feature (rf) trials?
  opts = list(n_folds = 5, pt_trials = 0, rf_trials = 0,
              paired = T, X_label = val_ID)

  return_vals2 <- cross_validation_mPLSDA(as.matrix(val_X), val_y,
                                          method, opts, n_trials = 100)

  val_plt <- visualize_validate(return_vals2)

  print(val_plt)

  #Average cross-validation scores across trials
  u2 <- rowMeans(sapply(return_vals2, unlist))[1]
  print(paste0(fig, " 5-fold cross validation mean score"))
  print(u2)
  u3 <- rowMeans(sapply(return_vals2,unlist))
  print(paste0(fig, " null models cross validation mean scores"))
  print(u3)
}
