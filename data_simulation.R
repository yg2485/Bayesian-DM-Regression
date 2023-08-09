# First, simulated dataset X with dim N*R
N <- 25
R <- 4
P <- 30
# Columns of X: age(1-70), race(0,1,2,3),gender(0,1) 0 as male, 
# phenotype(0,1) 0 as not sick

x_age <- rdunif(N, 10, 70)
x_race <- rbern(N, 0.7)
x_gender <- rbern(N, 0.3)
x_phenotype <- rbern(N, 0.8)
X <- cbind(range01(x_age), x_race, x_gender, x_phenotype)

Delta <- DATA$Delta
Beta <- DATA$Beta
Beta_0 <- DATA$Beta_0

setwd("~/Desktop/Research/DM Model/Research/Data from SW")
X_raw <- read.csv('Xmatrix.csv')
X = X_raw[,!(names(X_raw) %in% c("X"))]
X <- data.matrix(X, rownames.force = NA)
N <- dim(X)[1]
R <- dim(X)[2]
P <- 30

library(simcausal)
Delta <- matrix(rbern(R*P, 0.2), nrow = R, ncol = P)
Beta <- matrix(0, nrow = R, ncol = P)

for (i in 1:R){
  for (j in 1:P){
    if (Delta[i,j] == 1){
      p <- rbinom(1,1,0.5)
      beta <- runif(1, 1, 5)
      if (p == 1){
        beta <- -beta
      }
      Beta[i,j] <- beta
    }
  }
}

Beta_0 <- rtruncnorm(P, mean = 0, sd = 5, -3, 3)
B <- rbind(Beta_0, Beta)

X_design <- cbind(1, X)
A <- X_design %*% B
A <- exp(A)

# Simulated dataset Y with dim N*P
Y <- matrix(NA, nrow = N, ncol = P)
library(DirichletReg)
library(purrr)
for (i in 1:N){
  alpha <- A[i,]
  omega <- rdirichlet(1, alpha)
  total <- rdunif(1, 300, 7000)
  y <- rmultinom(1, total, omega) # Assume the total count for each sample is 1000
  Y[i,] <- y
}




T <-20000
d = 731
N <- dim(X)[1]
R <- dim(X)[2]
P <- dim(Y)[2]

result <- DMLN(X,Y,T = T, tau = 1, tau_beta_1 = 1, tau_beta_2 = 0.5, d = d) # run the main model

Beta_0_store <- result$Beta_0_store
Delta_store <- result$Delta_store
Beta_store <- result$Beta_store
acc_1 <- result$accept_1
acc_2 <- result$accept_2

deltaPPI <- matrix(NA, nrow = R, ncol = P)
for (i in 1:R){
  for (j in 1:P){
    deltaPPI[i,j] <- mean(Delta_store[which(Delta_store[,2]==i & Delta_store[,3]==j),4][(T/2):T])
  }
}


BetaPPI <- matrix(NA, nrow = R, ncol = P)
for (i in 1:R){
  for (j in 1:P){
    BetaPPI[i,j] <- mean(Beta_store[which(Beta_store[,2]==i & Beta_store[,3]==j),4][(T/2):T])
  }
}


# ========================================================================================
# ========================================================================================
# |||||||||||||||||||||||Simulated Data Analysis|||||||||||||||||||||||||||
# ========================================================================================
# ========================================================================================
cor(as.vector(BetaPPI), as.vector(Beta))

# Draw the ROC curve
delta_ROC(as.vector(deltaPPI), as.vector(abs(Delta)))

# Confusion Matrix
delta_pred <- matrix(NA, nrow = R, ncol = P)
for (i in 1:R){
  for (j in 1:P){
    if (deltaPPI[i,j] < 0.5){
      delta_pred[i,j] <- 0
    }else{
      delta_pred[i,j] <- 1
    }
  }
}

confusionMatrix(factor(as.vector(delta_pred)), factor(as.vector(abs(Delta))), mode = "everything", positive="1")

# For simulated data: correlation coefficient of real and posterior beta
Beta_0_post <- colMeans(Beta_0_store[(T/2):T,])

cor(Beta_0_post, Beta_0)



delta_ROC(as.vector(Beta_lm_c), as.vector(abs(Delta)))

#### Plot ROC curves in the same graph 
roc1 <- roc(as.vector(Delta), as.vector(deltaPPI))
roc2 <- roc(as.vector(Delta), as.vector(Beta_corr))
roc3 <- roc(as.vector(Delta), as.vector(Beta_lm_c))
plot(roc1, col = 'coral', lty = 2,main = 'ROC curves of multiple methods')
plot(roc2, col = "royalblue2", lty = 2, add = TRUE)
plot(roc3, col = 'darkgreen', lty = 2, add = TRUE)
legend(0.3, 0.2, legend = c('Bayesian PPI', 'Correlation Test', 'Linear Regression'), col = c('red', 'blue', 'green'),lty = 2, cex = 0.8)
legend(0.2, 0.45, title = 'AUC', legend = c('0.995', '0.785', '0.787'), col = c('red', 'blue', 'green'),lty = 2, cex = 0.8)

