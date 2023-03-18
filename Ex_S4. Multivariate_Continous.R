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
  list(x =X_E, model = "FIXED"),
  list(x =X_g, model = "BRR"),
  list(x =X_ge, model = "BRR")
)
head(Pheno)
# Retrieve the two continuos response variables
y <- SKM::to_matrix(Pheno[, c("GY", "PHR")])

# Folds generation for random cross validation
# Folds generation for random line validation
folds <- cv_random_line(lines =as.character(Pheno$Line),folds_number=5, testing_proportion =0.2)
#

# List for storing the individual predictions at each fold
Predictions <- list(GY = data.frame(), PHR = data.frame())

for (i in seq_along(folds)) {
  fold <- folds[[i]]
  
  # Set testing indices to NA in the response variables
  y_na <- y
  y_na[fold$testing, ] <- NA
  
  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::bayesian_model(
    X,
    y_na,
    
    iterations_number = 10000,
    burn_in = 5000
  )
  
  # Predict over testing set
  predictions <- predict(model)
  
  for (trait in names(Predictions)) {
    # Save the predictions along with environment and line information
    Predictions[[trait]] <- rbind(
      Predictions[[trait]],
      data.frame(
        Fold = i,
        Line = Pheno$Line[fold$testing],
        Env = Pheno$Env[fold$testing],
        Observed = y[fold$testing, trait],
        Predicted = predictions[[trait]]$predicted
      )
    )
  }
}

# Errors metrics for prediction performance
GY_Summaries <- SKM::gs_summaries(Predictions$GY)

# Print summaries by line, environment and fold
GY_Summaries$line
GY_Summaries$env
GY_Summaries$fold

# Errors metrics for prediction performance
PHR_Summaries <- SKM::gs_summaries(Predictions$PHR)

# Print summaries by line, environment and fold
PHR_Summaries$line
PHR_Summaries$env
PHR_Summaries$fold
