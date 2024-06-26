# A Regularized Bayesian Dirichlet-multinomial Regression Model for Integrating Single-cell-level Omics and Patient-level Clinical Study Data

The abundance of various cell types can vary significantly among patients with varying phenotypes and even those with the same phenotype. Recent scientific advancements provide mounting evidence that other clinical variables, such as age, gender, and lifestyle habits, can also influence the abundance of certain cell types. However, current methods for integrating single-cell-level omics data with clinical variables are inadequate. In this study, we propose a regularized Bayesian Dirichlet-multinomial regression framework to investigate the relationship between single-cell RNA sequencing data and patient-level clinical data. Additionally, the model employs a novel hierarchical tree structure to identify such relationships at different cell-type levels. Our model successfully uncovers significant associations between specific cell types and clinical variables across three distinct diseases: pulmonary fibrosis, COVID-19, and non-small cell lung cancer. This integrative analysis provides biological insights and could potentially inform clinical interventions for various diseases.

## Introduction
This is a repo for the paper "A Regularized Bayesian Dirichlet-multinomial Regression Model for Integrating Single-cell-level Omics and Patient-level Clinical Study Data". The repo contains the code for the simulation study and the real data analysis. 
### Required packages to run the code
The model was developed and tested under R 4.2.2
The following packages are required to run the code:
```r
library(purrr)
library(Rlab)
library(MCMCprecision)
library(ggplot2)
library(pROC)
library(reshape2)
library(caret)
library(dplyr)
library(Rcpp)
library(ape)
```

and the main function is in the file "DMLN.cpp".
```r
sourceCpp("./DMLN.cpp")
source("./helper_function.R")
```

## Real data analysis
Since the data is already preprocessed, we can directly import the data and do a basic data cleaning.


```r
X_raw <- read.csv("~/data/data2/X.csv", check.names=FALSE)
Y_raw <- read.csv("~/data/data2/Y.csv", check.names=FALSE)

X <- data.matrix(X_raw, rownames.force = NA)
Y <- data.matrix(Y_raw, rownames.force = NA)

# the following preprocessing steps is for the Covid-19 dataset (data2)
X <- X[,-c(4)]
ind <- which(rowSums(Y) > 200)
X <- X[ind,]
Y <- Y[ind,]

```

then we can run the DM model
```r
T <-20000
d = 731 # random seed
N <- dim(X)[1]
R <- dim(X)[2]
P <- dim(Y)[2]

result <- DMLN(X,Y,T = T, tau = 1, tau_beta_1 = 1, tau_beta_2 = 0.5, d = d) # run the main model

```
and save the results for further analysis
```r
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
```

then we can plot the results

### heatmap of the PPI of $\Delta$ and posterior $B$

```r
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

```
<p float="left">
    <img src='figure/data 2/data2-1.png' width='45%' height='45%' /> 
    <img src='figure/data 2/data2-2.png' width='45%' height='45%' /> 
</p>

### 95% credible interval of $B$
```r
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

```
this is the plot of the 95% credible interval of Beta

<img src='figure/data 2/2-1.png' width='50%' height='50%'> 
