#This will return Fig3A-B as well as cross-validation scores (may take a while to run)

#plotting options

opts_plot1 <- list(df_features = df_id,
                   loading_alpha = 1, # transparency for the loadings
                   score_alpha = 1, # transparency for the scores
                   color_features = "phenotype",
                   colors = my_colors,
                   y_name = "phenotype")

opts_plot3 <- list(df_features = df_features,
                   loading_alpha = 1, # transparency for the loadings
                   score_alpha = 1, # transparency for the scores
                   LV_ind = c(1,2), # which LVs to plot
                   n_LV = 2,
                   colors = my_colors,
                   y_name = "phenotype")

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

figure <- c("3A","3B")

for (fig in figure){
  if (fig == "3A"){
    #Z-score only serum data
    X_down <- data_corrected[which(fluid == "Serum"),]
    X_down <- as.matrix(scale(X_down, center = TRUE, scale = TRUE))
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    sep <- phenotype[which(fluid == "Serum")]
    sep <- droplevels(sep)
  } else if (fig == "3B"){
    #Z-score only joint fluid data
    X_down <- data_corrected[which(fluid == "Joint Fluid"),]
    X_down <- as.matrix(scale(X_down, center = TRUE, scale = TRUE))
    #X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    sep <- phenotype[which(fluid == "Joint Fluid")]
    sep <- droplevels(sep)
  }

  if (fig %in% c("3B")){
    #Take features selected in 80% of 100 rounds of LASSO
    opts_sel <- list(n_trials = 100, threshold = 0.8, return_count = FALSE)
  } else {
    #use Elastic Net feature selection & select an appropriate alpha with a grid search
    #because LASSO is too stringent for the serum model & doesn't work- while
    #Elastic Net also doesn't work as in it can't choose features to separate the groups, it still
    #works "better" than LASSO & is a better comparison if we're trying to compare
    #the best models possible for each context
    model <- caret::train(x = data.matrix(X_down), y = sep,method="glmnet",
                   trControl = trainControl("LOOCV"),tuneLength = 10)

    alpha_tuned <- model$bestTune$alpha

    opts_sel <- list(n_trials = 100, threshold = 0.7,
                     alpha = alpha_tuned, return_count = FALSE)
  }

  sel_features <- select_repeat(X_down, sep,
                                selector = select_lasso,
                                options = opts_sel)

  X_sel <- X_down[, sel_features]


  # Perform a PLS-DA using the selected features and plot the scores and loadings
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


for (fig in figure){
  if (fig == "3A"){
    #Z-score only serum data
    X_down <- data_corrected[which(fluid == "Serum"),]
    X_down <- as.matrix(scale(X_down, center = TRUE, scale = TRUE))
    X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    sep <- phenotype[which(fluid == "Serum")]
    sep <- droplevels(sep)
  } else if (fig == "BC"){
    #Z-score only joint fluid data
    X_down <- data_corrected[which(fluid == "Joint Fluid"),]
    X_down <- as.matrix(scale(X_down, center = TRUE, scale = TRUE))
    X_down <- X_down[,-which(colnames(X_down) %in% remove_ind)]
    sep <- phenotype[which(fluid == "Joint Fluid")]
    sep <- droplevels(sep)
  }

  if (fig %in% c("3B")){
    #Take features selected in 80% of 100 rounds of LASSO
    opts_sel <- list(n_trials = 100, threshold = 0.8, return_count = FALSE)
  } else {
    model <- train(x = data.matrix(X_down), y = sep, method="glmnet",
                   trControl = trainControl("LOOCV"),tuneLength = 10)

    alpha_tuned <- model$bestTune$alpha

    opts_sel <- list(n_trials = 100, threshold = 0.7,
                     alpha = alpha_tuned, return_count = FALSE)
  }

  select <- function(X, y) { return(select_repeat(X, y,
                                                  selector = select_lasso,
                                                  options = opts_sel)) }

  method = list(select = select,
                train = train_ropls,
                predict = predict_ropls,
                score = score_accuracy)

  #How many total trials, permutation (pt) trials, and random feature (rf) trials?
  opts = list(n_folds = 5, pt_trials = 100, rf_trials = 100,
              paired = F)

  return_vals2 <- cross_validation_unpaired(as.matrix(X_down), sep,
                                            method, opts, n_trials = 10)
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
