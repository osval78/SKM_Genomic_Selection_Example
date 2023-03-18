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


# Retrieve the two mixed response variables

Pheno$Cat=ifelse(Pheno$GY<quantile(Pheno$GY,probs=0.3),"Bottom",ifelse(Pheno$GY>quantile(Pheno$GY,probs=0.7),"Top","Middle"))
Pheno$Cat <- as.factor(Pheno$Cat)

y <- Pheno[, c("GY", "Cat")]

# Folds generation for stratified random cross validation
# it keeps the distribution of environment at each fold
folds <- SKM::cv_random_strata(
  data = Pheno$Env,
  folds_number = 5,
  testing_proportion = 0.25
)

# List for storing the individual predictions at each fold
Predictions <- list(DTHD = data.frame(), GY = data.frame())

for (i in seq_along(folds)) {
  fold <- folds[[i]]

  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training, ]
  y_testing <- y[fold$testing, ]

  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::deep_learning(
    X_training,
    y_training,

    # Tunable hyperparameters
    learning_rate = list(min = 0.001, max = 0.1),
    epochs_number = 50,
    batch_size = 32,
    layers = list(
      list(
        neurons_number = list(min = 20, max = 50),
        activation = "sigmoid",
        dropout = list(min = 0.2, max = 0.4)
      ),
      list(
        neurons_number = list(min = 20, max = 50),
        activation = "tanh",
        dropout = list(min = 0.2, max = 0.4)
      )
    ),

    # Tune configuration parameters
    tune_type = "Bayesian_optimization",
    tune_cv_type = "Random",
    tune_testing_proportion = 0.2,
    tune_folds_number = 2,
    tune_bayes_samples_number = 5,
    tune_bayes_iterations_number = 5,

    # Other algorithm's parameters
    optimizer = "adam",
    shuffle = TRUE,
    early_stop = TRUE,
    early_stop_patience = 10
  )

  # Predict over testing set
  predictions <- predict(model, X_testing)

  # Save the predictions along with environment and line information
  Predictions$GY <- rbind(
    Predictions$GY,
    data.frame(
      Fold = i,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing$GY,
      Predicted = predictions$GY$predicted
    )
  )

  # Save the predictions along with environment and line information
  Predictions$Cat <- rbind(
    Predictions$Cat,
    data.frame(
      Fold = i,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing$Cat,
      Predicted = predictions$Cat$predicted,
      predictions$Cat$probabilities
    )
  )
}

# Errors metrics for prediction performance
GY_Summaries <- SKM::gs_summaries(Predictions$GY)

# Print summaries by line, environment and fold
GY_Summaries$line
GY_Summaries$env
GY_Summaries$fold

# Errors metrics for prediction performance
Cat_Summaries <- SKM::gs_summaries(Predictions$Cat)

# Print summaries by line, environment and fold
Cat_Summaries$line
Cat_Summaries$env
Cat_Summaries$fold
