#####混合メンバシップ確率的ブロックモデル#####
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
k <- 10   #混合数
d <- 3000   #ノード数

#ノード構造を設定
theta <- extraDistr::rdirichlet(d, rep(0.2, 2*k))
omega <- rbeta(d, 2.5, 5.0)
gamma <- rbeta(2*k, 1.2, 1.5)
index_list <- list()
id_list <- list()

for(i in 1:d){
  #ノード確率を設定
  par <- matrix(theta[i, ], nrow=d-1, ncol=2*k) * theta[-i, ]
  z_rate <- par / rowSums(par)
  z_vec <- as.numeric(rmnom(d-1, 1, z_rate) %*% 1:(2*k))
  
  #ノードを生成
  index <- (1:d)[-i]
  flag <- rbinom(d-1, 1, omega[i] * gamma[z_vec])
  id_list[[i]] <- rep(i, d-1)
  index_list[[i]] <- index * flag
}

#リストを変換
node_vec0 <- unlist(index_list)
index_node <- which(node_vec0 > 0)
node_vec <- node_vec0[index_node]
id_vec <- unlist(id_list)[index_node]

#インデックスを作成
id_list1 <- list()
d_vec1 <- list()
id_list2 <- list()
d_vec2 <- list()
w <- rep(0, d)
for(i in 1:d){
  id_list1[[i]] <- which(id_vec==i)
  id_list2[[i]] <- which(node_vec==i)
  d_vec1[[i]] <- rep(1, length(id_list1[[i]]))
  d_vec2[[i]] <- rep(1, length(id_list2[[i]]))
  w[i] <- length(id_list1[[i]])
}
f <- sum(w)


##モデルに基づきデータを生成
#パラメータの事前分布
k0 <- k/2
alpha1 <- rep(0.15, k)
beta01 <- 1.0
beta02 <- 0.5
beta03 <- 0.5
beta04 <- 1.75

#パラメータを生成
phi <- matrix(0, nrow=k, ncol=k)
for(j in 1:k){
  index_phi <- cbind(rep(j, k), 1:k)
  for(l in 1:k){
    #ベータ分布からリンク確率を生成
    index <- index_phi[l, ]
    if(index[1]==index[2]){
      phi[index[1], index[2]] <- rbeta(1, beta01, beta02)
    } else {
      phi1 <- rbeta(1, beta03, beta04)     
      phi2 <- rbeta(1, (phi1)*k0, (1-phi1)*k0)
      phi[index[1], index[2]] <- phi1; phi[index[2], index[1]] <- phi2
    }
  }
}
phit <- phi
theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha1)
theta2 <- thetat2 <- extraDistr::rdirichlet(d, theta1 + 0.05)

#ノードごとにリンクを生成
y <- rep(0, f)
Z1_list <- list()
Z2_list <- list()

for(i in 1:d){
  #トピックを生成
  z1 <- rmnom(w[i], 1, theta1[i, ])
  z2 <- rmnom(w[i], 1, theta2[node_vec[id_list1[[i]]], ])
  z1_vec <- as.numeric(z1 %*% 1:k)
  z2_vec <- as.numeric(z2 %*% 1:k)
  
  #リンクを生成
  phi_vec <- rowSums(phi[z1_vec, ] * z2)
  y[id_list1[[i]]] <- rbinom(w[i], 1, phi_vec)
  
  #データを格納
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
}

#リストを変換
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
index_link <- which(y==1)




####マルコフ連鎖モンテカルロ法でMMSBを推定####
##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##事前分布の設定
alpha1 <- k
beta1 <- 0.1
beta2 <- 0.1

##パラメータの真値
theta1 <- thetat1
theta2 <- thetat2
phi <- phit

##パラメータの初期値
theta1 <- extraDistr::rdirichlet(d, rep(5.0, k))
theta2 <- extraDistr::rdirichlet(d, rep(5.0, k))
phi <- matrix(rbeta(k*k, 1.0, 2.5), nrow=k, ncol=k)
diag(phi) <- 0.5
Zi1 <- rmnom(f, 1, theta1[id_vec, ])
Zi2 <- rmnom(f, 1, theta2[node_vec, ]) 
z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)

##パラメータの格納用配列
THETA1 <- array(0, dim=c(d, k, R/keep))
THETA2 <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(d, k, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k)
SEG2 <- matrix(0, nrow=f, ncol=k)


##インデックスを作成
index10 <- id_vec[-index_link]
index11 <- id_vec[index_link]
index20 <- node_vec[-index_link]
index21 <- node_vec[index_link]

##対数尤度の基準値
#真値での対数尤度
LLbest <- sum(dbinom(y, 1, rowSums(phit[as.numeric(Z1 %*% 1:k), ] * Z2), log=TRUE))

#1パラメータモデルでの対数尤度
phi_mu <- mean(y)
LLst <- sum(dbinom(y, 1, phi_mu, log=TRUE))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##ノードのトピックを生成
  #トピックごとの期待尤度を計算
  Li1 <- matrix(0, nrow=f, ncol=k)
  Li2 <- matrix(0, nrow=f, ncol=k)
  
  for(j in 1:k){
    Li1[-index_link, j] <- theta1[index10, j] * (1-phi)[j, ][z2_vec[index20]]
    Li1[index_link, j] <- theta1[index11, j] * phi[j, ][z2_vec[index21]]
    Li2[-index_link, j] <- theta2[index20, j] * (1-phi)[, j][z1_vec[index10]]
    Li2[index_link, j] <- theta2[index21, j] * phi[, j][z1_vec[index11]]
  }
  
  #多項分布からノードのトピックを生成
  z1_rate <- Li1 / rowSums(Li1); z2_rate <-  Li2 / rowSums(Li2)   #潜在変数z
  Zi1 <- rmnom(f, 1, z1_rate); Zi2 <- rmnom(f, 1, z2_rate)   #トピックをサンプリング
  z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)
  Zi1_T <- t(Zi1); Zi2_T <- t(Zi2)

  
  ##パラメータをサンプリング
  #トピック分布をサンプリング
  wsum01 <- wsum02 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    wsum01[i, ] <- Zi1_T[, id_list1[[i]]] %*% d_vec1[[i]] 
    wsum02[i, ] <- Zi2_T[, id_list2[[i]]] %*% d_vec2[[i]]
  }
  wsum1 <- wsum01 + alpha1; wsum2 <- wsum02 + alpha1
  theta1 <- extraDistr::rdirichlet(d, wsum1)   #ディリクレ分布からthetaをサンプリング
  theta2 <- extraDistr::rdirichlet(d, wsum2)
  
  #リンク確率をサンプリング
  vsum0 <- t((1-y) * Zi1) %*% ((1-y) * Zi2) + beta1
  vsum1 <- t(y * Zi1) %*% (y * Zi2) + beta2
  phi <- matrix(rbeta(k*k, vsum1, vsum0), nrow=k, ncol=k)
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    PHI[, , mkeep] <- phi
  }
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG2 <- SEG2 + Zi2
  }
  
  if(rp%%disp==0){
    #対数尤度を計算
    phi_vec <- rep(0, f)
    z1_vec <- as.numeric(Zi1 %*% 1:k)
    z2_vec <- as.numeric(Zi2 %*% 1:k)
    for(i in 1:f){
      phi_vec[i] <- phi[z1_vec[i], z2_vec[i]]
    }
    #サンプリング結果の表示
    print(rp)
    print(c(sum(dbinom(y, 1, phi_vec, log=TRUE)), LLbest, LLst))
    print(round(cbind(phi, phit), 3))
  }
}


####サンプリング結果の可視化と要約####
#トピック分布のサンプリング結果
matplot(t(THETA1[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA1[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[10, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[15, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

#リンク確率のパラメータ
matplot(t(PHI[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(PHI[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")


