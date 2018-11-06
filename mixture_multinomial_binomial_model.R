#####混合多項-二項分布モデル#####
library(MASS)
library(gtools)
library(reshape2)
library(plyr)
library(dplyr)
library(knitr)
library(ggplot2)


####データの発生####
##モデルの設定
k <- 5   #セグメント数
n <- 500   #セグメントのサンプル数
N <- n * k   #総サンプル数
seg <- rep(1:k, rep(n, k))   #セグメント割当
ns1 <- rpois(N, 30)   #多項分布の個人ごとの頻度
ns2 <- rpois(N, 40)   #二項分布の個人ごとの頻度
th <- 10   #セグメント別のパラメータ数


##確率ベクトルを定義して、説明変数を発生
P1 <- matrix(0, nrow=k, ncol=th) 
P2 <- rep(0, k)
X1.list <- list()
X2 <- c()

for(j in 1:k){
  r <- ((j-1)*n+1):((j-1)*n+n)
  
  #確率を計算
  p <- rgamma(th, 0.8, 0.2)
  P1[j, ] <- p / sum(p)
  P2[j] <- runif(1, 0.15, 0.7)
  
  #応答変数を発生
  X1.list[[j]] <- t(apply(cbind(ns1[r], 1), 1, function(x) rmultinom(1, x[1], P1[j, ])))
  X2 <- c(X2, apply(cbind(ns2[r], 1), 1, function(x) rbinom(1, x[1], P2[j])))
}

X1 <- do.call(rbind, X1.list)


##データを結合して結果を集計
#多項データの要約
by(X1, seg, function(x) round(colMeans(x), 3))
by(X1, seg, function(x) summary(x))
by(X1, seg, function(x) round(colSums(x)/sum(x), 3))

#二項データの要約
by(X2, seg, function(x) round(mean(x), 3))
by(X2, seg, function(x) summary(x))
by(cbind(X2, ns2), seg, function(x) round(sum(x[1])/sum(x[2]), 3))


####EMアルゴリズムで混合多項-二項分布モデルを推定####
##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta1, theta2, n1, n2, X1, X2, r, k){
  LLind <- matrix(0, nrow=nrow(X1), ncol=k)
  for(i in 1:k){
    Li1 <- apply(cbind(n1, X1), 1, function(x) dmultinom(x[-1], x[1], theta1[i, ]))   #多項分布の尤度
    Li2 <- apply(cbind(n2, X2), 1, function(x) dbinom(x[-1], x[1], theta2[i]))   #二項分布の尤度
    Li <- Li1 * Li2
    LLind[, i] <- Li
  }
  
  LLho <- matrix(r, nrow=nrow(X1), ncol=k, byrow=T) * LLind   #観測データの尤度
  z <- LLho / matrix(apply(LLho, 1, sum), nrow=nrow(X1), ncol=k)   #潜在変数zの計算
  LLosum <- sum(log(apply(matrix(r, nrow=nrow(X1), ncol=k, byrow=T) * LLind, 1, sum)))   #観測データの対数尤度の計算
  rval <- list(LLob=LLosum, z=z, LL=LLind, Li1=Li1, Li2=Li2)
  return(rval)
}

#初期値の設定
iter <- 0
k <- 5   #セグメント数

##thetaの初期値の設定
theta1 <- matrix(0, nrow=k, ncol=th)
theta2 <- c()

for(j in 1:k){
  p <- runif(th, 0.1, 1.5)
  theta1[j, ] <- p / sum(p)
  theta2 <- c(theta2, runif(1, 0.25, 0.7))
}

##混合率rの初期値
r <- c(0.15, 0.15, 0.25, 0.2, 0.25)

#対数尤度の初期化
L <- LLobz(theta1=theta1, theta2=theta2, n1=ns1, n2=ns2, X1=X1, X2=X2, r=r, k=k)
LL1 <- L$LLob
z <- L$z
round(z, 3)

#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 1  

##EMアルゴリズム
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #Eステップの計算
  z <- L$z   #潜在変数zの出力
  
  #Mステップの計算と最適化
  #thetaの推定
  theta <- matrix(0, nrow=k, ncol=th)
  for(j in 1:k){
    #完全データの対数尤度からthetaの推定量を計算
    thetaseg1 <- colSums(matrix(z[, j], nrow=nrow(X1), ncol=th) * X1) / sum(z[, j] * ns1)   #重み付き多項分布の最尤推定
    thetaseg2 <- sum(z[, j]*X2) / sum(z[, j]*ns2)   #重み付き二項分布の最尤推定
    
    theta1[j, ] <-thetaseg1
    theta2[j] <- thetaseg2  
  }
  
  #混合率を推定
  r <- apply(z, 2, sum) / nrow(X1)
  
  #観測データの対数尤度を計算
  L <- LLobz(theta1=theta1, theta2=theta2, n1=ns1, n2=ns2, X1=X1, X2=X2, r=r, k=k)
  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}


####混合多項-二項分布モデルの推定結果####
round(theta1, 3)   #theta1の推定量
round(P1, 3)   #theta1の真の値
round(theta2, 3)   #theta1の推定量
round(P2, 3)   #theta1の真の値
round(r, 3)   #混合率の推定量
round(z, 3)   #個人別のセグメントへの所属確率

L$LLob   #観測データの対数尤度
-2*(L$LLob) + 2*(k*nrow(theta1)+length(theta2) + 1)   #AIC




