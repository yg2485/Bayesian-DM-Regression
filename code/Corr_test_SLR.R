library(reshape2)
######===========================================================######
################## Follow the steps to run the code ##################
######===========================================================######
# The following steps are for the correlation test and SLR for real data
# 1, import data. If you are interested in data 1, then only run the first part of "import data section"
# 2, run the correlation test and SLR. For data 1, run "Linear Regression for dataset 1", and
# if for data2 and 3, run "Linear Regression for dataset 2 and 3 (1 df lost in X)"
######===========================================================######




#################### Import data section ####################
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

#################### End of mport data section ####################





#################### Correlation Analysis section #####################

N <- dim(X)[1]
R <- dim(X)[2]
P <- dim(Y)[2]

Y_prop <- Y
for (i in 1:N){
  y_sum <- sum(Y[i,])
  for (j in 1:P){
    Y_prop[i,j] <- Y[i,j]/y_sum
  }
}


# Correlation Test
P_corr <- matrix(NA, nrow = R, ncol = P)
E_corr <- matrix(NA, nrow = R, ncol = P)
for (j in 1:P){
  for (i in 1:R){
    corr_test = cor.test(Y_prop[,j], X[,i])
    P_corr[i,j] <- corr_test$p.value
    E_corr[i,j] <- corr_test$estimate
  }
}

rownames(P_corr) = colnames(X)
colnames(P_corr) = colnames(Y)

data1 <- melt(P_corr)
data1_E <- as.vector(E_corr)
data1$value <- p.adjust(data1$value, method = 'fdr')
for (i in 1:(R*P)){
  if (!is.na(data1$value[i])){
    if(data1$value[i] > 0.05){
      data1_E[i] <- NA
    }
  }
}
data1$value2 <- data1_E
# delta_ROC(as.vector(Beta_corr), as.vector(abs(Delta)))



# Linear Regression for dataset 1
Beta_lm <- matrix(NA, nrow = R+1, ncol = P)
P_lm <- matrix(NA, nrow = R+1, ncol = P)
for (j in 1:P){
  y <- Y_prop[,j]
  data <- data.frame(y,X)
  l <- lm(y ~ ., data = data)
  summary_l <- summary(l)
  P_lm[,j] <- summary_l$coefficients[,4]
  Beta_lm[,j] <- summary_l$coefficients[,1]
}

nROW <- nrow(P_lm)
P_lm_c <- P_lm[2:nROW,]
Beta_lm_c <- Beta_lm[2:nROW,]

rownames(P_lm_c) = colnames(X)[1:(nROW-1)]
colnames(P_lm_c) = colnames(Y)

data2 <- melt(P_lm_c)
data2_E <- as.vector(Beta_lm_c)
data2$value <- p.adjust(data2$value, method = 'fdr')
for (i in 1:((nROW-1)*P)){
  if (!is.na(data2$value[i])){
    if(data2$value[i] > 0.05){
      data2_E[i] <- NA
    }
  }else{
    data2_E[i] <- NA
  }
}
data2$value2 <- data2_E



# # Linear Regression for dataset 2 and 3 (1 df lost in X)
# Beta_lm <- matrix(NA, nrow = R, ncol = P)
# P_lm <- matrix(NA, nrow = R, ncol = P)
# for (j in 1:P){
#   y <- Y_prop[,j]
#   data <- data.frame(y,X)
#   l <- lm(y ~ ., data = data)
#   summary_l <- summary(l)
#   P_lm[,j] <- summary_l$coefficients[,4]
#   Beta_lm[,j] <- summary_l$coefficients[,1]
# }
# 
# nROW <- nrow(P_lm)
# P_lm_c <- P_lm[2:nROW,]
# Beta_lm_c <- Beta_lm[2:nROW,]
# 
# rownames(P_lm_c) = colnames(X)[1:(nROW-1)]
# colnames(P_lm_c) = colnames(Y)
# 
# data2 <- melt(P_lm_c)
# data2_E <- as.vector(Beta_lm_c)
# data2$value <- p.adjust(data2$value, method = 'fdr')
# for (i in 1:((nROW-1)*P)){
#   if (!is.na(data2$value[i])){
#     if(data2$value[i] > 0.05){
#       data2_E[i] <- NA
#     }
#   }else{
#     data2_E[i] <- NA
#   }
# }
# data2$value2 <- data2_E

#################### End of Correlation Analysis section #####################




#################### Plot section ####################

ggplot(data1, aes(Var1, Var2, fill= value2)) + 
  geom_tile() +
  geom_text(aes(fill = value2, label = round(value2, 2)), size = 3) +
  scale_fill_gradient2(low="royalblue2", mid = "white",high="indianred1",na.value = "white")+
  labs(x = "Covariate", y = "Cell types",
       title ="Correlation coefficients with p-value <= 0.05")+
  theme(panel.background = element_blank(),
        legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text.x=element_text(angle = -60, hjust = 0, size = 20), 
        axis.text.y=element_text(hjust = 1, size = 13), 
        plot.title=element_text(hjust=0.5, size = 18))

ggplot(data2, aes(Var1, Var2, fill= value2)) + 
  geom_tile() +
  geom_text(aes(fill = value2, label = round(value2, 2)), size = 3) +
  scale_fill_gradient2(low="royalblue2", mid = "white",high="indianred1",na.value = "white")+
  labs(title ="Regression coefficients with p-value <= 0.05")+
  labs(x = "Covariate", y = "Cell types")+
  theme(panel.background = element_blank(),
        legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text.x=element_text(angle = -60, hjust = 0, size = 20), 
        axis.text.y=element_text(hjust = 1, size = 13), 
        plot.title=element_text(hjust=0.5, size = 18))

