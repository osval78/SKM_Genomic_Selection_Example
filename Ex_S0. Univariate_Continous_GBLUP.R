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

# Put the matrices in the expected format with the model
# to use with each one
X <- list(
  list(x =K_e, model = "BGBLUP"),
  list(x =K_g, model = "BGBLUP"),
  list(x =K_ge, model = "BGBLUP")
)
head(Pheno)
# Retrieve the two continuos response variables
y <- Pheno$GY
y
# Folds generation for random cross validation
# Folds generation for random line validation
folds <- cv_random_line(lines =as.character(Pheno$Line),folds_number=5, testing_proportion =0.2)
#

# List for storing the individual predictions at each fold
Predictions <-data.frame()

for (i in seq_along(folds)) {
#i=1 
 fold <- folds[[i]]
  
  # Set testing indices to NA in the response variables
  y_na <- y
  y_na[fold$testing] <- NA
  y_testing=y[fold$testing]
  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::bayesian_model(
    X,
    y_na,
    
    iterations_number = 10000,
    burn_in = 5000
  )
  
  # Predict over testing set
  predictions <- predict(model)
  
  # Save the predictions along with environment and line information
  Predictions <- rbind(
    Predictions,
    data.frame(
      Fold = i,
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
#Summaries$line
Summaries$env
#Summaries$fold
