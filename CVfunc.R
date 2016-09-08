# runCV is a helper function for running a machine learning 10-fold cross-validation predicting one ephys property using some list of metadata
#
# Inputs:
#   df - data.frame, NeuroElectro dataframe (NA's in the metadata have to be dealt with), ephys property will be filtered by this function
#   model - method to apply, Example: glmnet or randomForest
#   formula - String, Example: "rmp ~ NeuronName"
#   model_args - additional arguments to be passed on to the model (besides data and formula)
#   predict_args - additional arguments to be passed on to the fitting of model
#   splitNN - **important** expands NeuronName column into a matrix with each NeuronName being a column, required for randomForest that use NeuronName in model
#   unique - force each fold to have unique articles
#
# Output:
#   a data.frame with R^2 and Mean Squared Errors for each fold

get_R2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2 )
}

runCV <- function(df, model, formula, model_args, predict_args, splitNN = F, unique = F) {
  ep = gsub("\\s*~\\s*.*", "", formula)
  run_data = df[!is.na(df[, ep]),]
  run_data = run_data[sample(nrow(run_data)),]
  
  # Randomly assign folds to data rows
  run_data$fold = 0
  if (unique) {
    # Enforce unique Pmid's for each fold
    pmids <- unique(run_data$Pmid)
    for (i in 1:10) {
      split <- pmids[seq(i, length(pmids), 10)]
      run_data[run_data$Pmid %in% split,]$fold = i
    }
  } else {
    run_data$fold = rep(1:10, length.out = nrow(run_data))
  }
  
  if (splitNN & grepl("NeuronName", formula)) {
    rf_formula = gsub("NeuronName", "1", formula)
    run_data$NeuronName = droplevels(run_data$NeuronName)

    for(level in unique(run_data$NeuronName)){
      run_data[paste("NeuronName", make.names(level), sep = "_")] <- ifelse(run_data$NeuronName == level, 1, 0)
      rf_formula = paste(rf_formula, paste("NeuronName", make.names(level), sep = "_"), sep = "+")
    }
    run_data$NeuronName = NULL
    formula = gsub("1\\+", "", rf_formula)
  }
  
  result_r2 = c()
  result_mse = c()
  
  for (k in 1:10) {
    print(paste0("Fold #", k))
    
    fit = do.call(model, c(list(formula = as.formula(formula), data = run_data[run_data$fold != k,]), model_args) )
    pred = as.vector(do.call(predict, c(list(fit, newdata = run_data[run_data$fold == k,]), predict_args)))
    
    result_r2 = c(result_r2, get_R2(run_data[run_data$fold == k, ep], pred))
    result_mse = c(result_mse, mse(run_data[run_data$fold == k, ep], pred))
  }
  
  return(data.frame(fold = 1:10, r2 = result_r2, mse = result_mse))
}