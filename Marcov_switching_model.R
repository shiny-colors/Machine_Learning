#####マルコフ切り替えモデル#####
library(MASS)
library(MSwM) 
library(reshape2)
library(gtools)
library(dplyr)
library(ggplot2)
library(lattice)

####データの発生####
n <- 3000   #サンプル数
k1 <- 3   #切り替え数
k2 <- 2   #観測確率のパラメータ数


##初期確率の定義
Pf <- c(0.4, 0.3, 0.3)

##推移行列の定義
pr1 <- c(0.1, 0.7, 0.2)
pr2 <- c(0.2, 0.1, 0.7)
pr3 <- c(0.7, 0.2, 0.1)
Pr <- rbind(pr1, pr2, pr3)

##観測確率の定義
P <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  alpha <- runif(1, 0.4, 1)
  P[i, ] <- rdirichlet(1, rep(alpha, k2))
}

P <-rbind(c(0.9, 0.1), c(0.6, 0.4), c(0.1, 0.9))

##応答変数の発生
Z <- matrix(0, nrow=n, ncol=k1)
Y <- matrix(0, nrow=n, ncol=k2)

#潜在変数の初期値
Z[1, ] <- rmultinom(1, 1, Pf)
Y[1, ] <- rmultinom(1, 1, P[which.max(Z[1, ]), ])

#2回目以降の応答変数を逐次的に発生させる
for(i in 2:n){
  Z[i, ] <- rmultinom(1, 1, Pr[which.max(Z[i-1, ]), ])
  Y[i, ] <- rmultinom(1, 1, P[which.max(Z[i, ]), ])
}


####EMアルゴリズムでマルコフ切り替えモデルを推定####
##初期値の設定
#初期確率の設定
rho <- rep(0.25, k1)   

#マルコフ推移行列の初期値の設定
A <- matrix(0, nrow=k1, ncol=k1)
for(i in 1:k1){
  p_rand <- runif(k1, 0.1, 1)
  A[i, ] <- p_rand / sum(p_rand)
}

#観測モデルのパラメータ
B <- matrix(0, nrow=k1, ncol=k2)
for(i in 1:k1){
  p_rand <- runif(k2, 0.1, 1)
  B[i, ] <- p_rand / sum(p_rand)
}
y <- Y %*% 1:k2

#データの準備
y_delta <- array(0, dim=c(n, k1, k2))
for(i in 1:k2){
  y_delta[, , i] <- matrix(Y[, i], nrow=n, ncol=k1)
}

#パラメータの格納用配列
alpha <- matrix(0, nrow=n, ncol=k1)
alpha_s <- matrix(0, nrow=n, ncol=k1)
alpha_mu <- rep(0, n)
B_vec <- matrix(0, nrow=n, ncol=k1)
beta <- matrix(0, nrow=n, ncol=k1)
beta_s <- matrix(0, nrow=n, ncol=k1)

LL1 <- sum(log(apply(Y, 1, function(x) dmultinom(x, 1, B[1, ]))))   #対数尤度の初期値
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.01

L <- c()
for(i in 1:n){
  L <- c(L, dmultinom(Y[i, ], 1, P[which.max(Z[i, ]), ]))
}
LT <- sum(log(L))
LLl <- c()

####EMアルゴリズム(バウムウェルチアルゴリズム)で隠れマルコフモデルを推定####
while(abs(dl) >= tol){
  
  ##前向きアルゴリズムでalphaを推定
  B_vec[1, ] <- B[, y[1]]
  alpha[1, ] <- rho * B_vec[1, ]
  alpha_mu[1] <- 1 / sum(alpha[1, ])
  alpha_s[1, ] <- alpha_mu[1] * alpha[1, ] 
  
  for(i in 2:n){
    B_vec[i, ] <- B[, y[i]]
    alpha[i, ] <- alpha_s[i-1, ] %*% A * B_vec[i, ]
    alpha_mu[i] <- 1 / sum(alpha[i, ])
    alpha_s[i, ] <- alpha[i, ] * alpha_mu[i]
  }
  
  ##後ろ向きアルゴリズムでbetaを推定
  beta[n, ] <- 1
  beta_s[n, ] <- alpha_mu[n]
  
  for(i in n:2){
    beta[i-1, ] <- A %*% (B_vec[i, ] * beta_s[i, ])
    beta_s[i-1, ] <- beta[i-1, ] * alpha_mu[i-1] 
  }
  
  ##パラメータを更新
  #推移確率Aを更新
  for(i in 1:k1){
    A_vec <- matrix(A[i, ], nrow=n-1, ncol=k1, byrow=T)
    a11 <- matrix(alpha_s[1:(n-1), i], n-1, k1) * A_vec * B_vec[2:n, ] * beta_s[2:n, ]
    a12 <- matrix(alpha_s[1:(n-1), i], n-1, k1) * matrix(beta_s[1:(n-1), i], n-1, k1) / alpha_mu[1:(n-1)]
    a <- colSums(a11)/colSums(a12)
    A[i, ] <- a
  }
  
  
  #観測データのパラメータBを更新
  for(j in 1:k2){
    B[, j] <- colSums(y_delta[, , j] * alpha_s * beta_s / alpha_mu) / colSums(alpha_s * beta_s / alpha_mu)
  }
  
  #初期確率rhoを更新
  rho <- alpha_s[1, ] * beta_s[1, ] / alpha_mu[1]
  
  ##対数尤度を計算
  LL <- sum(log(1/alpha_mu))
  dl <- LL - LL1
  LL1 <- LL
  LLl <- c(LLl, LL)
  print(LL)
}

plot(1:length(LLl), LLl, type="l", main="対数尤度の変化", ylab="対数尤度", xlab="繰り返し数")
round(cbind(B, P), 3)
round(cbind(A, Pr), 3)
LT