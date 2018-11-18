#####mixed membership stochastic multiple block model#####
library(MASS)
library(bayesm)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(reshape2)
library(qrmtools)
library(slfm)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

set.seed(4965723)

####データの発生####
##データの設定
k <- 10   #グループ数
s <- 10   #リンクの種類数
hh <- 5000   #ユーザー数
item <- 3000   #アイテム数
company <- 500   #企業数
pt <- rpois(hh, rgamma(hh, 35, 0.25))   #ユーザーあたりの接触数
hhpt <- sum(pt)

##IDとインデックスの設定
#IDの設定
user_id <- rep(1:hh, pt)
pt_id <- as.numeric(unlist(tapply(1:hhpt, user_id, rank)))

#インデックスの設定
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(user_id==i)
}

##企業とアイテムの割当を生成
#セグメント割当を生成
topic <- 25
phi <- extraDistr::rdirichlet(topic, rep(0.5, item))
z <- as.numeric(rmnom(hh, 1,  extraDistr::rdirichlet(hh, rep(1.0, topic))) %*% 1:topic)

#多項分布からアイテムを生成
item_id_list <- list()
for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  item_id_list[[i]] <- as.numeric(rmnom(pt[i], 1, phi[z[user_id[user_list[[i]]]], ]) %*% 1:item)
}
item_id <- unlist(item_id_list)
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}
w <- unlist(lapply(item_list, length))

#アイテム割当からメーカーを生成
z <- rmnom(item, 1, as.numeric(extraDistr::rdirichlet(1, rep(1.5, company))))
index_z <- as.numeric(z[, colSums(z)!=0] %*% 1:sum(colSums(z)!=0))  
company_id <- left_join(data.frame(id=item_id), data.frame(id=1:item, no=index_z), by="id")$no   #idを設定
company <- length(unique(company_id))
company_list <- list()
for(j in 1:company){
  company_list[[j]] <- which(company_id==j)
}
company_id

#個別に和を取るためのスパース行列
user_dt <- sparseMatrix(user_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))
company_dt <- sparseMatrix(company_id, 1:hhpt, x=rep(1, hhpt), dims=c(company, hhpt))
item_dt <- sparseMatrix(item_id, 1:hhpt, x=rep(1, hhpt), dims=c(item, hhpt))


##パラメータの設定
#事前分布の設定
alpha11 <- rep(0.15, k)
alpha12 <- rep(0.1, k)
alpha13 <- rep(0.25, k)
alpha21 <- rep(0.2, s)
beta01 <- 5.0
beta02 <- 4.0
beta03 <- 1.25
beta04 <- 5.0

#パラメータを生成
omega <- matrix(0, nrow=k, ncol=k)
for(j in 1:k){
  index_omega <- cbind(rep(j, k), 1:k)
  for(l in 1:k){
    #ベータ分布からリンク確率を生成
    index <- index_omega[l, ]
    if(index[1]==index[2]){
      omega[index[1], index[2]] <- rbeta(1, beta01, beta02)
    } else {
      omega1 <- rbeta(1, beta03, beta04)     
      omega2 <- rbeta(1, (phi1)*(k/2), (1-phi1)*(k/2))
      omega[index[1], index[2]] <- omega1; omega[index[2], index[1]] <- omega2
    }
  }
}
omega[lower.tri(phi)] <- 0
omega <- phi + t(omega); diag(omega) <- diag(omega)/2



