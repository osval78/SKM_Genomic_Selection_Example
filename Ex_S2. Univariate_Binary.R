# Import SKM library
library(SKM)

# Load the dataset
load("Japonica_100.RData", verbose = TRUE)

# Oder Pheno by Env and Lines
Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
Markers=Markers
geno_order <- sort(rownames(Markers))
Markers=Markers[geno_order ,]
Markers_Scale=scale(Markers[,-1])
dim(Markers_Scale)
Markers_Scale=Markers_Scale[ , colSums(is.na(Markers_Scale))==0]
dim(Markers_Scale)
Geno=Markers_Scale%*%t(Markers_Scale)/ncol(Markers_Scale)

# Order Geno by Line
Geno <- Geno[geno_order, geno_order]
Markers_SVD=svd(Markers_Scale)
dim(Markers_SVD$u)
Markers_Red=Markers_SVD$u%*%diag(Markers_SVD$d)
dim(Markers_Red)
# Design matrices
Z_g<- model.matrix(~ 0 + Line, data = Pheno)
K_g=Z_g%*% Geno%*%t(Z_g) ###Linear kernel

# Env design matrix without the first column
X_E<- model.matrix(~ 0 + Env, data = Pheno)[, -1]
K_e=X_E%*%t(X_E)/ncol(X_E)

X_g <-Z_g%*%Markers_Red
K_ge <- K_e* K_g
# Interaction matrix
X_ge<- model.matrix(~ 0 + X_g:X_E, data = Pheno)
dim(X_ge)
# Bind all matrices in a single one
X <- cbind(X_E, X_g, X_ge)
# Retrieve the binary response variable
y=ifelse(Pheno$GY>quantile(Pheno$GY,probs=0.7),"Top","Bottom")
y <- as.factor(y)
# Folds generation for random line cross validation
folds <- SKM::cv_kfold_strata(
 y,
  k = 5
)

# Data frame for storing the individual predictions at each fold
Predictions <- data.frame()

for (i in seq_along(folds)) {
  fold <- folds[[i]]
  
  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training]
  y_testing <- y[fold$testing]
  
  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::support_vector_machine(
    X_training,
    y_training,
    
    kernel = "radial",
    gamma = list(min = 1/ncol(X), max = 0.1),
    cost = list(min = 0.5, max = 1),
    
    tune_type = "Bayesian_optimization",
    tune_bayes_samples_number = 5,
    tune_bayes_iterations_number = 5,
    
    tune_cv_type = "random",
    tune_folds_number = 5,
    tune_testing_proportion = 0.25,
    tune_loss_function = "roc_auc"
  )
  
  # Predict over testing set
  predictions <- predict(model, X_testing)
  
  # Save the predictions along with environment and line information
  Predictions <- rbind(
    Predictions,
    data.frame(
      Fold = i,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing,
      Predicted = predictions$predicted,
      predictions$probabilities
    )
  )
}

# Errors metrics for prediction performance
Summaries <- SKM::gs_summaries(Predictions)

# Print summaries by line, environment and fold
Summaries$line
Summaries$env
Summaries$fold
