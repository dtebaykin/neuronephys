# Created by: Dmitry Tebaykin
# Hosted on: https://github.com/dtebaykin/neuronephys
# Helper functions for running NeuroElectro prediction models

# This function returns R^2 value, inputs are 2 numerical vectors of equal length: observed values and predicted
get_R2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2 )
}

# Helper function for running a machine learning 10-fold cross-validation predicting one ephys property using some list of metadata
#
# Inputs:
#   df - data.frame, NeuroElectro dataframe (NA's in the metadata have to be dealt with), ephys property will be filtered by this function
#   model - method to apply, Example: glmnet or randomForest
#   formula - String, Example: "rmp ~ NeuronName"
#   model_args - additional arguments to be passed on to the model (besides data and formula)
#   predict_args - additional arguments to be passed on to the fitting of model
#   unique - force each fold to have unique articles
#   seed - enforce the splits to 'random' the same way each time, default - random seed between 1 and 10^5
#
# Output:
#   the original input data.frame with added columns: fold, pred
runCV <- function(df, model, formula, model_args = NULL, predict_args = NULL, unique = T, seed = sample(floor(runif(1, 1, 100001)), 1)) {
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
  
  # convert NeuronName column into a matrix where each NT is a column and rows are 1/0 based on occurrence in that article
  convertedList = convertNN(run_data, formula)
  
  run_data = convertedList[[1]]
  formula = convertedList[[2]]
  
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

# Convert NeuronName column into a matrix where each NeuronName is a column, adjust the formula accordingly
#
# Inputs:
#   df - data.frame with the NeuronName column
#   formula - Character, formula where NeuronName should be replaced with all the new columns being added
#
# Output:
#   list, first object - modified data.frame, second object - modified formula. Access via list[[1]] and list[[2]]
convertNN <- function(df, formula) {
  rf_formula = gsub("NeuronName", "1", formula)
  df$NeuronName = droplevels(df$NeuronName)
  
  for(level in unique(df$NeuronName)){
    df[paste("NeuronName", make.names(level), sep = "_")] <- ifelse(df$NeuronName == level, 1, 0)
    rf_formula = paste(rf_formula, paste("NeuronName", make.names(level), sep = "_"), sep = "+")
  }
  df$NeuronName = NULL
  formula = gsub("1\\+", "", rf_formula)
  
  return(list(df, formula))
}

# Run given model on the given dataset with the set parameters
# Inputs:
#   df - data.frame, input data to use for the model
#   model - Character, formula to use for the model
#   model_args - list, other parameters to be used
#
# Outputs:
#   model object
getModel <- function(df, model, formula, model_args = NULL) {
  df = df[!is.na(df[, ep]),]
  
  convertedList = convertNN(df, formula)
  run_data = convertedList[[1]]
  formula = convertedList[[2]]
  
  fit = do.call(model, c(list(formula = as.formula(formula), data = run_data), model_args) )
  return(fit)
}