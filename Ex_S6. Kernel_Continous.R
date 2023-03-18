###Example A6. Kernel model for an univariate continuous variable
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

X <- SKM::kernelize(
  X,
  kernel = "Sparse_Polynomial",
  rows_proportion = 0.75,
  gamma = 1 / ncol(X),
  coef0 = 0.5,
  degree = 4
)
# Retrieve the continuos response variable
y <- as.numeric(Pheno$GY)
str(y)

# Folds generation for k-fold cross validation
folds <- SKM::cv_kfold(nrow(Pheno), k = 5)

# Data frame for storing the individual predictions at each fold
Predictions <- data.frame()

for (fold_number in seq_along(folds)) {
  fold <- folds[[fold_number]]
  
  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training]
  y_testing <- y[fold$testing]
  
  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::random_forest(
    X_training,
    y_training,
    trees_number = c(100, 200),
    node_size = c(5, 2),
    
    tune_type = "Grid_search",
    tune_cv_type = "k_fold",
    tune_folds_number = 5,
    tune_loss_function = "mse"
  )
  
  # Predict over testing set
  predictions <- predict(model, X_testing)
  
  # Save the predictions along with environment and line information
  Predictions <- rbind(
    Predictions,
    data.frame(
      Fold = fold_number,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing,
      Predicted = predictions$predicted
    )
  )
}

# Errors metrics for prediction performance
Summaries <- SKM::gs_summaries(Predictions)

# Print summaries by line, environment and fold
Summaries$line
Summaries$env
Summaries$fold
