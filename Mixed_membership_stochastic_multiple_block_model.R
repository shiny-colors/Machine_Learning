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
k <- 15   #トピック数
s1 <- 10   #アイテムリンクの種類数
s2 <- 8   #企業リンクの種類数
hh <- 5000   #ユーザー数
item <- 3000   #アイテム数
company <- 1000   
pt1 <- rpois(hh, rgamma(hh, 25, 0.25))   #ユーザーあたりのアイテムへの接触数
pt2 <- rpois(hh, rgamma(hh, 30, 0.3))   #ユーザーあたりの企業への接触数
f1 <- sum(pt1)   #アイテムへの接触総数
f2 <- sum(pt2)   #企業への接触総数

##IDとインデックスの設定
#IDの設定
user_id1 <- rep(1:hh, pt1)
user_id2 <- rep(1:hh, pt2)
pt_id1 <- as.numeric(unlist(tapply(1:f1, user_id1, rank)))
pt_id2 <- as.numeric(unlist(tapply(1:f2, user_id2, rank)))

#インデックスの設定
user_list1 <- list()
user_list2 <- list()
for(i in 1:hh){
  user_list1[[i]] <- which(user_id1==i)
  user_list2[[i]] <- which(user_id2==i)
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
  item_id_list[[i]] <- as.numeric(rmnom(pt1[i], 1, phi[z[user_id1[user_list1[[i]]]], ]) %*% 1:item)
}
item_id <- unlist(item_id_list)
item_list <- list()
for(j in 1:item){
  item_list[[j]] <- which(item_id==j)
}
w <- unlist(lapply(item_list, length))

#個別に和を取るためのスパース行列
user_dt1 <- sparseMatrix(user_id1, 1:f1, x=rep(1, f1), dims=c(hh, f1))
user_dt2 <- sparseMatrix(user_id2, 1:f2, x=rep(1, f2), dims=c(hh, f2))
item_dt <- sparseMatrix(item_id, 1:f1, x=rep(1, f1), dims=c(item, f1))


##パラメータの設定
#事前分布の設定
alpha11 <- rep(0.15, k)
alpha12 <- rep(0.15, k)
alpha21 <- rep(0.05, company)
alpha22 <- rep(0.2, s1)
alpha23 <- rep(0.2, s2)
beta01 <- 5.0; beta02 <- 4.0
beta03 <- 1.25; beta04 <- 5.0


##すべてのデータが出現するまでデータの生成を続ける
rp <- 0
repeat { 
  rp <- rp + 1
  print(rp)
  
  #ディリクレ分布からトピック分布を生成
  theta1 <- extraDistr::rdirichlet(hh, alpha11)
  theta2 <- extraDistr::rdirichlet(item, alpha12)
  
  #パラメータを生成
  gamma1 <- extraDistr::rdirichlet(k, alpha22)
  gamma2 <- extraDistr::rdirichlet(k, alpha23)
  phi <- extraDistr::rdirichlet(k, alpha21)
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
        omega2 <- rbeta(1, (omega1)*(k/2), (1-omega1)*(k/2))
        omega[index[1], index[2]] <- omega1; omega[index[2], index[1]] <- omega2
      }
    }
  }
  omega[lower.tri(omega)] <- 0
  omega <- omega + t(omega); diag(omega) <- diag(omega)/2
  
  #出現確率が低いphiの要素を入れ替える
  index <- which(colMaxs(phi) < (k*5)/f2)
  for(j in 1:length(index)){
    phi[as.numeric(rmnom(1, 1, extraDistr::rdirichlet(1, alpha11)) %*% 1:k), index[j]] <- (k*5)/f2
  }

  ##応答変数を生成
  RX <- matrix(0, nrow=hh, ncol=company)
  Z11_list <- Z12_list <- Z21_list <- list()
  y_list <- r_list <- link1_list <- link2_list <- list()

  for(i in 1:hh){
    #多項分布からユーザーのトピックを生成
    Z11 <- rmnom(pt1[i], 1, theta1[i, ])
    Z12 <- rmnom(pt2[i], 1, theta1[i, ])
    z_vec11 <- as.numeric(Z11 %*% 1:k)
    z_vec12 <- as.numeric(Z12 %*% 1:k)
    
    #多項分布からアイテムのトピックを生成
    Z21 <- rmnom(pt1[i], 1, theta2[item_id[user_list1[[i]]], ])
    z_vec21 <- as.numeric(Z21 %*% 1:k)
    
    #トピックからリンクとリンクの種類を生成
    omega_vec <- as.numeric((omega[z_vec11, ] * Z21) %*% rep(1, k))
    y <- rbinom(pt1[i], 1, omega_vec)
    link1 <- rmnom(pt1[i], 1, gamma1[z_vec11, ])
    
    #トピックから企業とリンクの種類を生成
    r <- rmnom(pt2[i], 1, phi[z_vec12, ])
    link2 <- rmnom(pt2[i], 1, gamma2[z_vec12, ])
    
    #データを格納
    RX[i, ] <- colSums(r)
    Z11_list[[i]] <- Z11
    Z12_list[[i]] <- Z12
    Z21_list[[i]] <- Z21
    y_list[[i]] <- y
    r_list[[i]] <- r
    link1_list[[i]] <- link1
    link2_list[[i]] <- link2
  }
  #break条件
  if(min(colSums(RX)) > 0){
    break
  }
}





