#####Nested Latent bariable Gaussian Distribution Model#####
options(warn=2)
library(MASS)
library(mclust)
library(flexmix)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(mvtnorm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

####多変量正規分布の乱数を発生させる関数を定義####
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}


##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.01, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの生成####
##データの設定
k1 <- 5
k2 <- 7
d <- 1000   #データ数
w <- rpois(d, rgamma(d, 30, 0.7))
f <- sum(w)
v <- 5   #変数数

#IDの設定
d_id <- rep(1:d, w)
n_id <- c()
for(i in 1:d){
  n_id <- c(n_id, 1:w[i])
}

##パラメータの設定
#ディリクレ分布のパラメータ
alpha1 <- rep(0.4, k1)
alpha2 <- rep(0.5, k2)


#多変量正規分布のパラメータ
Cov <- Covt <- array(0, dim=c(v, v, k1))
Mu <- Mut <- matrix(0, nrow=k1, ncol=v)
for(j in 1:k1){
  corr <- corrM(v, -0.8, 0.9, 0.01, 0.2)
  Cov[, , j] <- Covt[, , j] <- covmatrix(v, corr, 2^2, 4^2)$covariance
  Mu[j, ] <- Mut[j, ] <- runif(v, 10.0, 30.0)
}

#固定パラメータ
beta0 <- betat0 <- mvrnorm(k2, rep(0, v), diag(9.0, v))

#ディリクレ分布のパラメータを生成
theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha1)
theta2 <- thetat2 <- extraDistr::rdirichlet(k1, alpha2)


##LHAモデルに基づきデータを生成
Data_list <- Z1_list <- Z2_list <- list()

for(i in 1:d){
  #基底を生成
  z1 <- rmnom(w[i], 1, theta1[i, ])
  z1_vec <- as.numeric(z1 %*% 1:k1)
  
  #固定パラメータの割当を生成
  z2 <- rmnom(w[i], 1, theta2[z1_vec, ])
  z2_vec <- as.numeric(z2 %*% 1:k2)
  
  #多変量正規分布から観測データを生成
  mu <- Mu[z1_vec, ] + beta0[z2_vec, ]   #平均パラメータ
  y <- matrix(0, nrow=w[i], ncol=v)
  for(j in 1:w[i]){
    y[j, ] <- mvrnorm(1, mu[j, ], Cov[, , z1_vec[j]])
  }
  
  #データを格納
  Data_list[[i]] <- y
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
}

#リストを変換
Data <- do.call(rbind, Data_list)
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)


####変分ベイズEMアルゴリズムでHLAを推定####
##多変量正規分布の尤度関数
dmv <- function(x, mean.vec, S, S_det, S_inv){
  LLo <- (2*pi)^(-nrow(S)/2) * S_det^(-1/2) *
    exp(-1/2 * (x - mean.vec) %*% S_inv %*% (x - mean.vec))
  return(LLo)
}

##アルゴリズムの設定
R <- 2000
keep <- 2  
iter <- 0
burnin <- 200/keep
disp <- 10

##事前分布の設定
#多変量正規分布の事前分布
mu0 <- rep(0, v)
sigma0 <- 100
sigma0_inv <- 1/sigma0

#逆ウィシャート分布の事前分布
nu <- v + 1
V <- nu * diag(v)
inv_V <- solve(V)

#ディリクレ分布の事前分布
alpha <- 1

##パラメータの真値
theta1 <- thetat1
theta2 <- thetat2
mu <- Mu
Cov <- Covt
beta0 <- betat0


##初期値の設定
theta1 <- extraDistr::rdirichlet(d, rep(10.0, k1))
theta2 <- extraDistr::rdirichlet(k1, rep(10.0, k2))
mu <- mvrnorm(k1, rep(mean(Data), v), diag(2^2, v))
Cov <- array(diag(2^2, v), dim=c(v, v, k1))
beta0 <- mvrnorm(k2, rep(0, v), diag(2^2, v))

##パラメータの格納用配列
THETA1 <- array(0, dim=c(d, k1, R/keep))
THETA2 <- array(0, dim=c(k1, k2, R/keep))
MU <- array(0, dim=c(k1, v, R/keep))
BETA <- array(0, dim=c(k2, v, R/keep))
COV <- array(0, dim=c(v, v, k1, R/keep))
SEG <- matrix(0, nrow=f, ncol=k1*k2)


#インデックスを作成
index_k1 <- matrix(1:(k1*k2), nrow=k1, ncol=k2, byrow=T)
index_column <- rep(1:k1, rep(k2, k1))
d_list <- d_vec <- list()
for(i in 1:d){
  d_list[[i]] <- which(d_id==i)
  d_vec[[i]] <- rep(1, length(d_list[[i]]))
}


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##潜在変数zをサンプリング
  #多変量正規分布の対数尤度を計算
  Li <- matrix(0, nrow=f, ncol=k1*k2)
  for(i in 1:k1){
    for(j in 1:k2){
      Li[, index_k1[i, ][j]] <- dmvnorm(Data, mu[i, ] + beta0[j, ] , Cov[, , i], log=TRUE)
    }
  }
  
  #潜在変数zの事前分布の設定
  theta1_z <- log(theta1)[d_id, index_column]
  theta2_z <- matrix(as.numeric(t(log(theta2))), nrow=f, ncol=k1*k2, byrow=T)
  
  
  #潜在変数zの計算
  LLi <- theta1_z + theta2_z + Li   #潜在変数zの対数尤度
  z_par <- exp(LLi - rowMaxs(LLi))   #尤度に変換
  z_rate <- z_par / rowSums(z_par)   #潜在変数z
  
  #多項分布から潜在変数をサンプリング
  Zi <- rmnom(f, 1, z_rate)
  z_vec <- as.numeric(Zi %*% 1:(k1*k2))
  
  #パターンごとに潜在変数を割当
  Zi1 <- matrix(0, nrow=f, ncol=k1)
  Zi2 <- matrix(0, nrow=f, ncol=k2)
  for(j in 1:k1) {Zi1[, j] <- rowSums(Zi[, index_k1[j, ]])}
  for(j in 1:k2) {Zi2[, j] <- rowSums(Zi[, index_k1[, j]])}
  z1_vec <- as.numeric(Zi1 %*% 1:k1)
  z2_vec <- as.numeric(Zi2 %*% 1:k2)
  
  #混合率を更新
  r <- colSums(Zi) / f
  
  
  ##多変量正規分布のパラメータと固定パラメータを更新
  #平均ベクトルを更新
  index1 <- list()
  for(j in 1:k1){
    index1[[j]] <- which(Zi1[, j]==1)
    mu_par <- colSums(Data[index1[[j]], ] - beta0[as.numeric(Zi[index1[[j]], index_k1[j, ]] %*% 1:k2), ])
    mu_mean <- mu_par / (length(index1[[j]]) + sigma0_inv)   #平均パラメータ
    mu_cov <- Cov[, , j] / (1+length(index1[[j]]))   #平均ベクトルの分散共分散行列
    mu[j, ] <- mvrnorm(1, mu_mean, mu_cov)   #多変量正規分布より平均ベクトルをサンプリング
  }
  
  #固定パラメータを更新
  for(j in 1:k2){
    weighted_cov <- matrix(0, nrow=v, ncol=v)
    index2 <- which(Zi2[, j]==1)
    mu_par <- colSums(Data[index2, ] - mu[as.numeric(Zi[index2, index_k1[, j]] %*% 1:k1), ])
    mu_mean <- mu_par / (length(index2) + sigma0_inv)   #平均パラメータ
    for(l in 1:k1){
      weighted_cov <-  weighted_cov <- Cov[, , l] * mean(Zi[index2, index_k1[, j]][, l])
    }
    mu_cov <- weighted_cov / (1+length(index2))   #固定パラメータの分散共分散行列
    beta0[j, ] <- mvrnorm(1, mu_mean, mu_cov)   #多変量正規分布より固定パラメータをサンプリング
  }
  
  #分散共分散行列の変分事後平均を推定
  for(j in 1:k1){
    Vn <- nu + length(index1[[j]])
    er <- Data[index1[[j]], ] - mu[z1_vec[index1[[j]]], ] - beta0[as.numeric(Zi[index1[[j]], index_k1[j, ]] %*% 1:k2), ]
    R_par <- solve(V) + t(er) %*% er
    Cov[, , j] <- rwishart(Vn, solve(R_par))$IW   #逆ウィシャート分布から分散共分散行列をサンプリング
  }
  
  ##潜在変数の割当確率を更新
  #基底の分布を更新
  Zi1_T <- t(Zi1)   #潜在変数zの転置行列
  wsum0 <- matrix(0, nrow=d, ncol=k1)
  for(i in 1:d){
    wsum0[i, ] <- Zi1_T[, d_list[[i]]] %*% d_vec[[i]]
  }
  wsum1 <- wsum0 + alpha   #ディリクレ分布のパラメータ
  theta1 <- extraDistr::rdirichlet(d, wsum1)   #ディリクレ分布からサンプリング
  
  
  #固定パラメータの割当分布の更新
  for(j in 1:k1){
    wsum2 <- colSums(Zi[, index_k1[j, ]]) + alpha   #ディリクレ分布のパラメータ
    theta2[j, ] <- extraDistr::rdirichlet(1, wsum2)   #ディリクレ分布からサンプリング
  }
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    MU[, , mkeep] <- mu
    BETA[, , mkeep] <- beta0
    COV[, , , mkeep] <- Cov
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      SEG <- SEG + Zi
    }
    
    if(rp%%disp==0){
      #サンプリング結果を確認
      print(rp)
      print(sum(log(rowSums(exp(LLi)))))
      print(round(cbind(mu, Mut), 3))
      print(round(cbind(beta0, betat0), 3))
    }
  }
}

####サンプリング結果の可視化と要約####
burnin <- 500/keep
RS <- R/keep

##サンプリング結果の可視化
#多変量正規分布のパラメータをプロット
matplot(t(MU[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(MU[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(MU[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(BETA[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

#基底分布のパラメータをプロット
matplot(t(THETA1[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA1[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA1[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA2[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

##パラメータの事後平均を推定
#多変量正規分布のパラメータの事後平均
round(t(apply(MU[, , burnin:RS], c(1, 2), mean)), 3)   #平均ベクトル
round(t(apply(BETA[, , burnin:RS], c(1, 2), mean)), 3)   #固定パラメータ
for(j in 1:k1){
  print(apply(COV[, , j, ], c(1, 2), mean))
}

#基底分布のパラメータの事後平均
round(apply(THETA1[, , burnin:RS], c(1, 2), mean), 3)
round(apply(THETA2[, , burnin:RS], c(1, 2), mean), 3)




