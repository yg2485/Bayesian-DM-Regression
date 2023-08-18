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
sourceCpp("~/Desktop/Research/DM Model/Research/DMLN.cpp")
source("~/Desktop/Research/DM Model/Research/helper_function.R")


####### Import data #######
#===============================================================================
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

X <- X[,-c(6)]
ind <- which(rowSums(Y) > 200)
X <- X[ind,]
Y <- Y[ind,]

################################################################################
#===============================================================================
################################################################################

N <- dim(X)[1]
R <- dim(X)[2]
P <- dim(Y)[2]


############################## dendrogram_object ##############################

# Data 1
load("/Users/yanghongguo/Desktop/Research/DM Model/Research/data/dendrogram/dendrogram_1.RData")
plot(dendrogram_1)

# Data 2
load("/Users/yanghongguo/Desktop/Research/DM Model/Research/data/dendrogram/dendrogram_2.RData")
plot(dendrogram_2)

# Data 3
load("/Users/yanghongguo/Desktop/Research/DM Model/Research/data/dendrogram/dendrogram_3.RData")
plot(dendrogram_3)


CC <- cophenetic(dendrogram_1)
CC <- cophenetic(dendrogram_2)
CC <- cophenetic(dendrogram_3) # get the dist matrix from dendrogram_object


# hierarchical clustering object
hc <- hclust(CC)
plot(hc)



# ct1 <- cutree(hc, k = 3)
label_all <- hc$labels
label_all <- gsub("_", ".", label_all)
label_all <- gsub("\\+", ".", label_all)
label_all <- gsub("-", ".", label_all)
label_all <- gsub(",", ".", label_all)
hc$labels <- label_all
# check identity
identical(sort(unique(label_all)), sort(colnames(Y)))
setdiff(unique(label_all), unique(colnames(Y)))



######### get the leaves of all the parent nodes #########
node_index <- list()

deal_m <- function(m, save_leaf, save_node){
  if (m[1] < 0 & m[2] < 0){
    save_leaf <- c(save_leaf, -m[1], -m[2])
  }else if(m[1] < 0 & m[2] > 0){
    save_leaf <- c(save_leaf, -m[1])
    save_node <- c(save_node, m[2])
  }else if (m[1] > 0 & m[2] < 0){
    save_leaf <- c(save_leaf, -m[2])
    save_node <- c(save_node, m[1])
  }else{
    save_node <- c(save_node, m[1], m[2])
  }
  return (list(save_leaf = save_leaf, save_node = save_node))
}


for (n in 1:(P-2)){
  
  m <- hc$merge[n,]
  res <- deal_m(m, c(), c())
  save_leaf <- res$save_leaf
  save_node <- res$save_node
  while (length(save_node) > 0){
    m <- hc$merge[save_node[1],]
    save_node <- save_node[-1]
    res <- res <- deal_m(m, save_leaf, save_node)
    save_leaf <- res$save_leaf
    save_node <- res$save_node
  }
  
  node_index[[n]] <- save_leaf
  
}
#########


# count how many times each node appear in all the layers
count_each_node <- vector("numeric", length = P-2) 


# generate Y_matrix for all layers
for (l in 1:(P-2)){
  # l represents the layer, from bottom to top
  # e.g. l=1 means only the first node merged
  num_merge <- l
  used <- c() # leaves that appears on one parent node
  name_node <- c()
  name_leaf <- c()
  Y_mat <- matrix(nrow = N, ncol = 0)
  for (s in l:1){
    label_num <- node_index[[s]]
    if (!all(label_num %in% used)){
      count_each_node[s] <- count_each_node[s] + 1
      label_sub <- label_all[label_num]
      Y_hat <- rowSums(Y[,label_sub])
      Y_mat <- cbind(Y_mat, Y_hat)
      used <- c(used, label_num)
      name_node <- c(name_node, paste0("node_", s))
    }
  }
  alll <- seq(1:P)
  left <- setdiff(alll, used)
  if (length(left) > 0 ){
    name_leaf <- paste("leave_", left, sep = "")
  }
  
  Y_mat <- cbind(Y_mat, Y[,label_all[left]])
  colnames(Y_mat) <-c(name_node, name_leaf)
  
  # define name of Y_mat as Y_mat_l
  obj_name <- paste("Y_mat_", l, sep = "")
  
  # Create the object and assign it a value
  assign(obj_name, Y_mat)
}


##### run the DMLN #####
T = 20000
d = 731
for (l in 1:(P-2)){
  message(sprintf("This is for layer %d", l))
  assign("Y_mat", get(paste0("Y_mat_", l)))
  Y_mat <- as.matrix(Y_mat)
  res <- DMLN(X, Y_mat, T = T, tau = 1, tau_beta_1 = 1, tau_beta_2 = 0.5, d = d)
  
  Beta_0_store <- res$Beta_0_store
  Delta_store <- res$Delta_store
  Beta_store <- res$Beta_store
  
  N <- dim(X)[1]
  R <- dim(X)[2]
  Q <- dim(Y_mat)[2]
  
  deltaPPI <- matrix(NA, nrow = R, ncol = Q)
  for (i in 1:R){
    for (j in 1:Q){
      deltaPPI[i,j] <- mean(Delta_store[which(Delta_store[,2]==i & Delta_store[,3]==j),4][(T/2):T])
    }
  }
  
  BetaPPI <- matrix(NA, nrow = R, ncol = Q)
  for (i in 1:R){
    for (j in 1:Q){
      BetaPPI[i,j] <- mean(Beta_store[which(Beta_store[,2]==i & Beta_store[,3]==j),4][(T/2):T])
    }
  }
  
  colnames(deltaPPI) <- colnames(Y_mat)
  colnames(BetaPPI) <- colnames(Y_mat)
  # define name of deltaPPI and BetaPPI
  obj_name_1 <- paste("deltaPPI_", l, sep = "")
  obj_name_2 <- paste("BetaPPI_", l, sep = "")
  
  # Create the object and assign it a value
  assign(obj_name_1, deltaPPI)
  assign(obj_name_2, BetaPPI)
  
  # define name of Beta_store
  obj_name_3 <- paste("Beta_store_", l, sep = "")
  assign(obj_name_3, Beta_store)
  
}

# Finally! Build the Delta and Beta matrix for each pair of parent node and covariate
deltaPPI_all <- matrix(NA, nrow = R, ncol = (P-2))
BetaPPI_all <- matrix(NA, nrow = R, ncol = (P-2))
Beta_025_all <- matrix(NA, nrow = R, ncol = (P-2))
Beta_975_all <- matrix(NA, nrow = R, ncol = (P-2))

all_nodes <- c()

for (l in 1:(P-2)){
  num_of_appear <- count_each_node[l]
  name_this_node <- paste0("node_", l)
  row_mean_delta <- matrix(NA, nrow = R, ncol = 0)
  row_mean_beta <- matrix(NA, nrow = R, ncol = 0)
  quan <- matrix(NA, nrow = 0, ncol = R)
  
  for (t in 1:num_of_appear){
    name_Y_mat <- paste("Y_mat_", l + t - 1, sep = "")
    name_deltaPPI <- paste("deltaPPI_", l + t - 1, sep = "")
    name_BetaPPI <- paste("BetaPPI_", l + t - 1, sep = "")
    name_Beta_store <- paste("Beta_store_", l + t - 1, sep = "")
    
    # assign deltaPPI_l, BetaPPI_l, Beta_store_l to temp1,2,3
    assign("temp1", get(name_deltaPPI))
    assign("temp2", get(name_BetaPPI))
    assign("temp3", get(name_Beta_store))
    assign("tempY", get(name_Y_mat))
    
    index_in_Y <- which(colnames(tempY) == name_this_node)
    
    # deltaPPI and BetaPPI
    row_mean_delta <- cbind(row_mean_delta, temp1[,name_this_node])
    row_mean_beta <- cbind(row_mean_beta, temp2[,name_this_node])
    
    # Quantiles
    temp_store <- matrix(NA, nrow = (T/2), ncol = R)
    for (rr in 1:R){
      temp_store[,rr] <- temp3[which(temp3[,2] == rr & temp3[,3] == index_in_Y),4][(T/2+1):T]
    }
    quan <- rbind(quan, temp_store)
  }
  row_mean_delta <- rowSums(row_mean_delta)/num_of_appear
  row_mean_beta <- rowSums(row_mean_beta)/num_of_appear
  deltaPPI_all[,l] <- row_mean_delta
  BetaPPI_all[,l] <- row_mean_beta
  
  Beta_025_all[,l] <- apply(quan, 2, quantile, probs = 0.025)
  Beta_975_all[,l] <- apply(quan, 2, quantile, probs = 0.975)
  
  # save all node names
  all_nodes <- c(all_nodes, name_this_node)
  # colnames(deltaPPI_all)[l] <- name_this_node
  # colnames(BetaPPI_all)[l] <- name_this_node
}


# assign rowname to deltaPPI_all and BetaPPI_all

colnames(deltaPPI_all) <- all_nodes
colnames(BetaPPI_all) <- all_nodes
rownames(deltaPPI_all) = colnames(X)
rownames(BetaPPI_all) = colnames(X)


#======================= Use ggplot to plot ============================
threshold <- BayFDR(as.vector(deltaPPI_all), 0.05)

BetaPPI2 <- BetaPPI_all
for (i in 1:R){
  for (j in 1:(P-2)){
    if(deltaPPI_all[i,j] < threshold){
      BetaPPI2[i,j] <- NA
    }
  }
}


data1 <- reshape2::melt(deltaPPI_all)
data2 <- reshape2::melt(BetaPPI2)

value2 <- array(NA, (R*(P-2)))
for (i in 1:(R*(P-2))){
  if (data1$value[i] >= threshold){
    value2[i] <- data1$value[i]
  }else{
    value2[i] <- NA
  }
}
data1$value2 <- value2


ggplot(data1, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  geom_text(aes(fill = value2, label = round(value2, 2)), size = 3, col="white") +
  scale_fill_gradient2(low="white", high="gray30", na.value = "white")+
  labs(title = expression(paste("Heatmap of the PPI of ", Delta)))+
  theme(panel.background = element_blank(),
        legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text.x=element_text(angle = -60, hjust = 0, size = 12), 
        axis.text.y=element_text(hjust = 0, size = 10), 
        plot.title=element_text(hjust=0.5))

ggplot(data2, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  geom_text(aes(fill = value, label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low="royalblue2", mid = "white",high="indianred1",na.value = "white")+
  labs(title = expression(paste("Heatmap of covariate effect on cell type (", beta, ")")))+
  theme(panel.background = element_blank(),
        legend.title=element_blank(),
        axis.title = element_blank(),
        axis.text.x=element_text(angle = -60, hjust = 0, size = 12), 
        axis.text.y=element_text(hjust = 0, size = 10), 
        plot.title=element_text(hjust=0.5))+
  guides(fill = guide_colourbar(ticks = FALSE))




#================================||||||||||||||||||||||||||==================================
#=============================== Plot 95% interval for Beta =================================
#================================||||||||||||||||||||||||||==================================

row_name <- rownames(deltaPPI_all)
col_name <- colnames(deltaPPI_all)

name_pair <- c()
quantiles_l <- c()
quantiles_u <- c()
means <- c()
counts <- 0
for (i in 1:R){
  for (j in 1:(P-2)){
    if(deltaPPI_all[i,j] >= threshold){
      new_name <- paste0(row_name[i],'__',col_name[j])
      name_pair <- c(name_pair, new_name)
      quantiles_l <- c(quantiles_l, Beta_025_all[i,j])
      quantiles_u <- c(quantiles_u, Beta_975_all[i,j])
      means <- c(means,BetaPPI_all[i,j])
    }
  }
}
post_Beta <- cbind(quantiles_l, quantiles_u, means)
rownames(post_Beta) <- name_pair
post_Beta <- data.frame(post_Beta)
post_Beta <- post_Beta[order(post_Beta$means,decreasing = FALSE),]
post_Beta$rowname <- as.character(rownames(post_Beta))
#Then turn it back into a factor with the levels in the correct order
post_Beta$rowname <- factor(post_Beta$rowname, levels=unique(post_Beta$rowname))
post_Beta$group <- ifelse(post_Beta$means > 0, "positive", "negative")

par(mar = c(2,2,2,2))
# new version
j <- ggplot(post_Beta)
j + geom_errorbar(aes(x = rowname, ymin = quantiles_l, ymax = quantiles_u,color = as.factor(group)), width = 0.2, size = 1)+
  scale_color_manual(values = c("positive" = "indianred1", "negative" = "royalblue")) +
  geom_point(aes(rowname, means), shape = 16, size = 2) +
  geom_hline(linetype = "dotted", size = 0.5, aes(yintercept = 0))+
  # labs(title = expression(paste("95% credible intervals of ", beta, " for tree")))+
  labs(x = "Cell Type and Covariate Pairs", y = expression(paste("Values of ", beta)), color = NULL)+
  theme( panel.background = element_blank(),
         axis.title = element_text(size = 15),
         axis.line = element_line(size = 0.2, color = "black"),
         axis.text.x=element_text(hjust = 0, size = 15), 
         axis.text.y=element_text(size = 15), 
         plot.title=element_text(hjust=0.5, size = 18),
         plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
         legend.position="none")+
  coord_flip(ylim = c(-2, 2))

