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
#   the original input data.frame with added columns: fold, pred

# This function returns R^2 value, inputs are 2 numerical vectors of equal length: observed values and predicted
get_R2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2 )
}

# Main function of the script - returns input data.frame with pred and fold columns
runCV <- function(df, model, formula, model_args, predict_args, splitNN = T, unique = T, seed = 42) {
  set.seed(seed)
  
  # Misc
  ep = gsub("\\s*~\\s*.*", "", formula)
  df$index = 1:nrow(df)
  
  # Filter out NA's in the ephys property we are predicting and reshuffle whatever is left
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
  
  # Convert NeuronName column into a matrix where each NeuronName is a column, adjust the formula accordingly
  rf_formula = gsub("NeuronName", "1", formula)
  run_data$NeuronName = droplevels(run_data$NeuronName)

  for(level in unique(run_data$NeuronName)){
    run_data[paste("NeuronName", make.names(level), sep = "_")] <- ifelse(run_data$NeuronName == level, 1, 0)
    rf_formula = paste(rf_formula, paste("NeuronName", make.names(level), sep = "_"), sep = "+")
  }
  run_data$NeuronName = NULL
  formula = gsub("1\\+", "", rf_formula)
  
  # Run the model for each of 10 folds
  run_data$pred = NA
  for (k in 1:10) {
    print(paste0("Fold #", k))
    
    fit = do.call(model, c(list(formula = as.formula(formula), data = run_data[run_data$fold != k,]), model_args) )
    run_data[run_data$fold == k,]$pred = as.vector(do.call(predict, c(list(fit, newdata = run_data[run_data$fold == k,]), predict_args)))
  }
  
  # Sort the predictions based on original data.frame
  run_data = run_data[order(run_data[,"index"]),]
  
  # Save folds
  df$fold = NA
  df[!is.na(df[, ep]),]$fold = run_data$fold
  
  # Save predictions
  df$pred = NA
  df[!is.na(df[, ep]),]$pred = run_data$pred
  
  # Delete index
  df$index = NULL
  
  return(df)
}