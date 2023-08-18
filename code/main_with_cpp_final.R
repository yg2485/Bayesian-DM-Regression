library(purrr)
library(Rlab)
library(MCMCprecision)
library(ggplot2)
library(pROC)
library(reshape2)
library(caret)
library(dplyr)
library(Rcpp)
sourceCpp("~/Desktop/Research/DM Model/Research/DMLN.cpp")
source("~/Desktop/Research/DM Model/Research/helper_function.R")

#################### Simulated data set ####################
# refer to "data_simulation.R"

#################### Real data set ####################

##### data 1 #####
X_raw <- read.csv("~/Desktop/Research/DM Model/Research/data/final/data1/X.csv", check.names=FALSE)
Y_raw <- read.csv("~/Desktop/Research/DM Model/Research/data/final/data1/Y.csv", check.names=FALSE)

X <- data.matrix(X_raw, rownames.force = NA)
Y <- data.matrix(Y_raw, rownames.force = NA)


##### data 2 #####
X_raw <- read.csv("~/Desktop/Research/DM Model/Research/data/final/data2/X.csv", check.names=FALSE)
Y_raw <- read.csv("~/Desktop/Research/DM Model/Research/data/final/data2/Y.csv", check.names=FALSE)

X <- data.matrix(X_raw, rownames.force = NA)
Y <- data.matrix(Y_raw, rownames.force = NA)

X <- X[,-c(4)]
ind <- which(rowSums(Y) > 200)
X <- X[ind,]
Y <- Y[ind,]


# colnames(X) <- c("age", "sex", "sc-seq.platform", "control", "convalescence&moderate",
#                  "convalescence&severe", "progression&moderate", "progression&severe")

##### data 3 #####
X_raw <- read.csv("~/Desktop/Research/DM Model/Research/data/final/data3/X.csv", check.names=FALSE)
Y_raw <- read.csv("~/Desktop/Research/DM Model/Research/data/final/data3/Y.csv", check.names=FALSE)

X <- data.matrix(X_raw, rownames.force = NA)
Y <- data.matrix(Y_raw, rownames.force = NA)
# colnames(X_raw) <- c("age", "sex", "ever.smoker", "tumor.stage.advanced", "tumor.stage.early", "tumor.stage.noncancer")

X <- X[,-c(6)]

ind <- which(rowSums(Y) > 200)
X <- X[ind,]
Y <- Y[ind,]

# ========================================================================================
# ========================================================================================
############################# Run MCMC Model #############################
# ========================================================================================
# ========================================================================================

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
# |||||||||||||||||||||||||||Real Data Analysis||||||||||||||||||||||||||||||||
# ========================================================================================
# ========================================================================================

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

######################### HEATMAP ###########################
#======================= Use ggplot to plot ============================
threshold <- BayFDR(as.vector(deltaPPI), 0.05)

BetaPPI2 <- BetaPPI
for (i in 1:R){
  for (j in 1:P){
    if(deltaPPI[i,j] < threshold){
      BetaPPI2[i,j] <- NA
    }
  }
}
rownames(BetaPPI2) = colnames(X)
colnames(BetaPPI2) = colnames(Y)

rownames(deltaPPI) = colnames(X)
colnames(deltaPPI) = colnames(Y) 

data1 <- reshape2::melt(deltaPPI)
data2 <- reshape2::melt(BetaPPI2)

value2 <- array(NA, (R*P))
for (i in 1:(R*P)){
  if (data1$value[i] >= threshold){
    value2[i] <- data1$value[i]
  }else{
    value2[i] <- NA
  }
}
data1$value2 <- value2


ggplot(data1, aes(Var1, Var2, fill= value)) + 
  # geom_tile(color = "gray",linetype = 2, size = .01) +
  geom_tile() +
  geom_text(aes(fill = value2, label = round(value2, 2)), size = 3, col="white") +
  scale_fill_gradient2(low="white", high="gray30", na.value = "white")+
  # labs(x = "Covariates", y = "Cell types",
  #      title = expression(paste("Heatmap of the PPI of ", Delta)))+
  # labs(title = expression(paste("Heatmap of the PPI of ", Delta)))+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x=element_text(angle = -60, hjust = 0, size = 18), 
        axis.text.y=element_text(hjust = 1, size = 13), 
        plot.title=element_text(hjust=0.5, size = 18))+
  guides(fill = guide_colorbar(title = NULL))

ggplot(data2, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  geom_text(aes(fill = value, label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low="royalblue2", mid = "white",high="indianred1",na.value = "white")+
  # labs(title = expression(paste("Heatmap of covariate effect on cell type (", beta, ")")))+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.x=element_text(angle = -60, hjust = 0, size = 18), 
        axis.text.y=element_text(hjust = 1, size = 13), 
        plot.title=element_text(hjust=0.5, size = 18))+
  guides(fill = guide_colourbar(title = NULL, ticks = FALSE))



#================================||||||||||||||||||||||||||==================================
#=============================== Plot 95% interval for Beta =================================
#================================||||||||||||||||||||||||||==================================

row_name <- colnames(X)
col_name <- colnames(Y)

name_pair <- c()
quantiles <- c()
means <- c()
counts <- 0
for (i in 1:R){
  for (j in 1:P){
    if(deltaPPI[i,j] >= threshold){
      counts <- counts + 1
      new_name <- paste0(row_name[i],'__',col_name[j])
      name_pair <- c(name_pair, new_name)
      sub_beta <- Beta_store[which(Beta_store[,2] == i & Beta_store[,3] == j),4]
      quantiles <- cbind(quantiles, quantile(sub_beta[(T/2):T], c(0.025, 0.975)))
      means <- c(means,BetaPPI[i,j])
    }
  }
}
post_Beta <- cbind(t(quantiles),means)
rownames(post_Beta) <- name_pair
post_Beta <- data.frame(post_Beta)
post_Beta <- post_Beta[order(post_Beta$means,decreasing = FALSE),]
post_Beta$rowname <- as.character(rownames(post_Beta))
#Then turn it back into a factor with the levels in the correct order
post_Beta$rowname <- factor(post_Beta$rowname, levels=unique(post_Beta$rowname))
post_Beta$group <- ifelse(post_Beta$means > 0, "positive", "negative")


# new version
# par(mar = c(3,3,3,3))
j <- ggplot(post_Beta)
j + geom_errorbar(aes(x = rowname, ymin = X2.5., ymax = X97.5.,color = as.factor(group)), width = 0.2, size = 1)+
  scale_color_manual(values = c("positive" = "indianred1", "negative" = "royalblue")) +
  geom_point(aes(rowname, means), shape = 16, size = 2) +
  geom_hline(linetype = "dotted", size = 0.5, aes(yintercept = 0))+
  # labs(title = expression(paste("95% credible intervals of ", beta)), color = NULL)+
  labs(x = "cell type and covariate pairs", y = expression(paste("values of ", beta)), color = NULL)+
  theme( panel.background = element_blank(),
         axis.title = element_text(size = 15),
         axis.line = element_line(size = 0.2, color = "black"),
         axis.text.x=element_text(hjust = 0, size = 15), 
         axis.text.y=element_text(size = 15), 
         plot.title=element_text(hjust=0.5, size = 18),
         plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
         legend.position="none")+
  coord_flip(ylim = c(-2, 2))



#================================||||||||||||||||||||||||||==================================
#=============================== Plot pie chart/ bar plot =================================
#================================||||||||||||||||||||||||||==================================
DATA <- readRDS(file = '~/Desktop/Research/DM Model/Research/data/final result/Real Data Result3.RData')
Beta_0_store <- DATA$Beta_0_store
Delta_store <- DATA$Delta_store
Beta_store <- DATA$Beta_store
X <- DATA$X
Y <- DATA$Y
N <- DATA$N
P <- DATA$P
R <- DATA$R
T <- DATA$T


deltaPPI <- matrix(NA, nrow = R, ncol = P)

for (i in 1:R){
  for (j in 1:P){
    deltaPPI[i,j] <- mean(Delta_store[which(Delta_store[,2]==i & Delta_store[,3]==j),4][(T/2):T])
  }
}

threshold <- BayFDR(as.vector(deltaPPI), 0.05)

BetaPPI <- matrix(NA, nrow = R, ncol = P)
for (i in 1:R){
  for (j in 1:P){
    if(deltaPPI[i,j] < threshold){
      BetaPPI[i,j] <- 0
    }else{
      BetaPPI[i,j] <- mean(Beta_store[which(Beta_store[,2]==i & Beta_store[,3]==j),4][(T/2):T])
    }
  }
}
Beta_0_post <- colMeans(Beta_0_store[(T/2):T,])

Beta_post <- rbind(Beta_0_post, BetaPPI)

labels <- colnames(Y)

#=========== now deal with each dataset ===========#

### dataset 1
# female
# for no disease (disease = 0)
x_n <- c()
x_n[1] <- 1
x_n[2] <- 0
x_n[3] <- mean(X[,"age"])
x_n[4] <- 1
x_n[5] <- 0

# for with disease (disease = 1)
x_y <- c()
x_y[1] <- 1
x_y[2] <- 0
x_y[3] <- mean(X[,"age"])
x_y[4] <- 1
x_y[5] <- 1


alpha_n <- exp(x_n %*% Beta_post)
alpha_y <- exp(x_y %*% Beta_post)

diff_1 <- c(alpha_n[5]/sum(alpha_n) ,alpha_y[5]/sum(alpha_y))
diff_2 <- c(alpha_n[8]/sum(alpha_n), alpha_y[8]/sum(alpha_y))

Y_norm <- t(apply(Y, 1, function(x) x/sum(x)))

# add category
category <- ifelse(X[,4] == 1, "With Disease", "No Disease")

Y_norm <- data.frame(Y_norm)
Y_norm$category = category

mid = barplot(diff_1, names.arg = c("No Disease", "With Disease"), main = "Epi (basal)",
        ylim = c(-0.005,0.1), col=c("seagreen3","salmon1"), border=F, las = 1, 
        cex.axis=1.5, cex.names=2, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)), 0.4), Y_norm[,5][order(Y_norm$category)], pch=19)


mid = barplot(diff_2, names.arg = c("No Disease", "With Disease"), main = "Epi (KRT5-/KRT17+)",
        ylim = c(-0.005,0.1), col=c("seagreen3","salmon1"), border=F, las = 1, 
        cex.axis=1.5, cex.names=2, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,8][order(Y_norm$category)], pch=19)







### dataset 2
# female  # platform 10*3'
x_1 <- rep(0, 8) # control group
x_2 <- rep(0, 8) # Convalescence.moderate
x_3 <- rep(0, 8) # Convalescence.severe
x_4 <- rep(0, 8) # Progression.moderate
x_5 <- rep(0, 8) # Progression.severe

x_1[1] = x_2[1] = x_3[1] = x_4[1] = x_5[1] = 1
x_1[2] = x_2[2] = x_3[2] = x_4[2] = x_5[2] = mean(X[,"age"])
x_2[5] = x_3[6] = x_4[7] = x_5[8] = 1

alpha_1 <- exp(x_1 %*% Beta_post)
alpha_2 <- exp(x_2 %*% Beta_post)
alpha_3 <- exp(x_3 %*% Beta_post)
alpha_4 <- exp(x_4 %*% Beta_post)
alpha_5 <- exp(x_5 %*% Beta_post)

# CD4
diff_1 <- c(alpha_1[2]/sum(alpha_1), alpha_2[2]/sum(alpha_2), alpha_3[2]/sum(alpha_3), alpha_4[2]/sum(alpha_4), alpha_5[2]/sum(alpha_5))

# CD8
diff_2 <- c(alpha_1[3]/sum(alpha_1), alpha_2[3]/sum(alpha_2), alpha_3[3]/sum(alpha_3), alpha_4[3]/sum(alpha_4), alpha_5[3]/sum(alpha_5))

# Mono
diff_3 <- c(alpha_1[9]/sum(alpha_1), alpha_2[9]/sum(alpha_2), alpha_3[9]/sum(alpha_3), alpha_4[9]/sum(alpha_4), alpha_5[9]/sum(alpha_5))

# NK
diff_4 <- c(alpha_1[11]/sum(alpha_1), alpha_2[11]/sum(alpha_2), alpha_3[11]/sum(alpha_3), alpha_4[11]/sum(alpha_4), alpha_5[11]/sum(alpha_5))

df <- data.frame(name = c("Control", "Convalescence.moderate", "Convalescence.severe", "Progression.moderate", "Progression.severe"),
                 diff_1 = diff_1,
                 diff_2 = diff_2,
                 diff_3 = diff_3,
                 diff_4 = diff_4)


Y_norm <- t(apply(Y, 1, function(x) x/sum(x)))

# add category
category <- rep("", nrow(Y_norm))
for (i in 1:nrow(Y_norm)){
  if (X[i,4] == 1){
    category[i] <- "Convalescence.moderate"
  }else if(X[i,5] == 1){
    category[i] <- "Convalescence.severe"
  }else if(X[i,6] == 1){
    category[i] <- "Progression.moderate"
  }else if(X[i,7] == 1){
    category[i] <- "Progression.severe"
  }else{
    category[i] <- "Control"
  }
}

Y_norm <- data.frame(Y_norm)
Y_norm$category = category

par(mar = c(15,3,2,2))
mid = barplot(df$diff_1, ylim = c(-0.01,0.7), names.arg = df$name, main = "CD4",
              col=c("seagreen3","salmon1","steelblue3","orchid2","lightskyblue2"), border=F,
              las = 2, cex.axis=1.5, cex.names=1.4, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)), 0.4), Y_norm$CD4[order(Y_norm$category)], pch=20)

mid = barplot(df$diff_2, ylim = c(-0.01,0.8), names.arg = df$name, main = "CD8",
              col=c("seagreen3","salmon1","steelblue3","orchid2","lightskyblue2"), border=F,
              las = 2, cex.axis=1.5, cex.names=1.4, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)), 0.4), Y_norm$CD8[order(Y_norm$category)], pch=20)

mid = barplot(df$diff_1, ylim = c(-0.01,1), names.arg = df$name, main = "Mono",
              col=c("seagreen3","salmon1","steelblue3","orchid2","lightskyblue2"), border=F,
              las = 2, cex.axis=1.5, cex.names=1.4, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)), 0.4), Y_norm$Mono[order(Y_norm$category)], pch=20)

mid = barplot(df$diff_1, ylim = c(-0.01,0.5), names.arg = df$name, main = "NK",
              col=c("seagreen3","salmon1","steelblue3","orchid2","lightskyblue2"), border=F,
              las = 2, cex.axis=1.5, cex.names=1.4, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)), 0.4), Y_norm$NK[order(Y_norm$category)], pch=20)








### dataset 3
# female
x_adv <- rep(0, 6)
x_early <- rep(0, 6)
x_non <- rep(0, 6)

x_adv[1] = x_early[1] = x_non[1] = 1
x_adv[2] = x_early[2] = x_non[2] = mean(X[,"age"])
x_adv[4] = x_early[4] = x_non[4] = 1 # smoked
x_adv[5] = x_early[6] = 1 

alpha_adv <- as.numeric(exp(x_adv %*% Beta_post))
alpha_early <- as.numeric(exp(x_early %*% Beta_post))
alpha_non <- as.numeric(exp(x_non %*% Beta_post))

# Macro (alv)
diff_1 <- c(alpha_adv[10]/sum(alpha_adv), alpha_early[10]/sum(alpha_early), alpha_non[10]/sum(alpha_non))

# T (CD4+)
diff_2 <- c(alpha_adv[12]/sum(alpha_adv), alpha_early[12]/sum(alpha_early), alpha_non[12]/sum(alpha_non))

# T (CD8+)
diff_3 <- c(alpha_adv[13]/sum(alpha_adv), alpha_early[13]/sum(alpha_early), alpha_non[13]/sum(alpha_non))

# Neutrophil
diff_4 <- c(alpha_adv[16]/sum(alpha_adv), alpha_early[16]/sum(alpha_early), alpha_non[16]/sum(alpha_non))

# T (reg)
diff_5 <- c(alpha_adv[19]/sum(alpha_adv), alpha_early[19]/sum(alpha_early), alpha_non[19]/sum(alpha_non))

# Mono (non)
diff_6 <- c(alpha_adv[20]/sum(alpha_adv), alpha_early[20]/sum(alpha_early), alpha_non[20]/sum(alpha_non))

df <- data.frame(name = c("Advanced-stage", "Early-stage", "Non-cancer"),
                 diff_1 = diff_1,
                 diff_2 = diff_2,
                 diff_3 = diff_3,
                 diff_4 = diff_4,
                 diff_5 = diff_5,
                 diff_6 = diff_6)

Y_norm <- t(apply(Y, 1, function(x) x/sum(x)))

# add category
category <- rep("", nrow(Y_norm))
for (i in 1:nrow(Y_norm)){
  if (X[i,4] == 1){
    category[i] <- "Advanced-stage"
  }else if(X[i,5] == 1){
    category[i] <- "Early-stage"
  }else{
    category[i] <- "Non-cancer"
  }
}

Y_norm <- data.frame(Y_norm)
Y_norm$category = category

par(mar = c(10,3.6,2,2))
mid = barplot(df$diff_1, ylim = c(-0.01,0.8), names.arg = df$name, main = "Macro (alv)",
              col=c("seagreen3","salmon1","steelblue3"), border=F,
              las = 1, cex.axis=1.5, cex.names=1.25, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,10][order(Y_norm$category)], pch=20)


mid = barplot(df$diff_2, ylim = c(-0.01,0.6), names.arg = df$name, main = "T (CD4+)",
              col=c("seagreen3","salmon1","steelblue3"), border=F,
              las = 1, cex.axis=1.5, cex.names=1.25, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,12][order(Y_norm$category)], pch=20)


mid = barplot(df$diff_3, ylim = c(-0.01,0.7), names.arg = df$name, main = "T (CD8+)",
              col=c("seagreen3","salmon1","steelblue3"), border=F,
              las = 1, cex.axis=1.5, cex.names=1.25, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,13][order(Y_norm$category)], pch=20)


mid = barplot(df$diff_4, ylim = c(-0.01,0.5), names.arg = df$name, main = "Neutrophil",
              col=c("seagreen3","salmon1","steelblue3"), border=F,
              las = 1, cex.axis=1.5, cex.names=1.25, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,16][order(Y_norm$category)], pch=20)


mid = barplot(df$diff_5, ylim = c(-0.01,0.3), names.arg = df$name, main = "T (reg)",
              col=c("seagreen3","salmon1","steelblue3"), border=F,
              las = 1, cex.axis=1.5, cex.names=1.25, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,19][order(Y_norm$category)], pch=20)


mid = barplot(df$diff_6, ylim = c(-0.01,0.85), names.arg = df$name, main = "Mono (non)",
              col=c("seagreen3","salmon1","steelblue3"), border=F,
              las = 1, cex.axis=1.5, cex.names=1.25, cex.main = 2)
points(jitter(rep(mid, table(Y_norm$category)),0.4), Y_norm[,20][order(Y_norm$category)], pch=20)

