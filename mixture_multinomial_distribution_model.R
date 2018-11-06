#####混合多項分布モデル(ログサムexp利用)#####
library(MASS)
library(vcd)
library(gtools)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(caret)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(lattice)

#set.seed(6439)

####データの発生####
##データの設定
hh <- 20000   #ユーザー数
k <- 10   #セグメント数
v <- 400   #アイテム
s <- rpois(hh, 20)   #購買数
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
f <- sum(s)

##パラメータの設定
theta <- thetat <- extraDistr::rdirichlet(1, rep(10, k))
phi <- phit <- extraDistr::rdirichlet(k, rep(0.3, v))

##多項分布よりデータを生成
Z_list <- list()
Data <- matrix(0, nrow=hh, ncol=v)

for(i in 1:hh){
  #セグメント割当を生成
  z <- rmnom(1, 1, theta)
  z_vec <- as.numeric(z %*% 1:k)
  Z_list[[i]] <- z
  
  #購買データを生成
  w <- rmnom(s[i], 1, phi[z_vec, ])
  Data[i, ] <- colSums(w)
}
Z <- do.call(rbind, Z_list)

####EMアルゴリズムで混合多項分布を推定####
const <- lfactorial(s) - rowSums(lfactorial(Data))   #多項分布の密度関数の対数尤度の定数
r <- colMeans(Z)

##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(Data, phi, r, const, hh, k){
  
  #多項分布の対数尤度
  log_phi <- log(t(phi))
  LLi <- const + Data %*% log_phi
  
  #logsumexpの尤度
  LLi_max <- matrix(apply(LLi, 1, max), nrow=hh, ncol=k)
  r_matrix <- matrix(r, nrow=hh, ncol=k, byrow=T)
  
  #割当確率のパラメータを設定
  expl <- r_matrix * exp(LLi - LLi_max)
  expl_log <- log(expl)
  expl_max <- matrix(log(max(expl[1, ])), nrow=hh, ncol=k)
  z <- exp(expl_log - (log(rowSums(exp(expl_log - expl_max))) + expl_max))   #セグメント割当確率
  
  #観測データの対数尤度
  r_log <- matrix(log(r), nrow=hh, ncol=k, byrow=T)
  LLosum <- sum(log(rowSums(exp(r_log + LLi))))   #観測データの対数尤度
  rval <- list(LLob=LLosum, z=z, LL=LLi)
  return(rval)
}

##パラメータの初期値]
#phiの初期値
alpha0 <- colSums(Data) / sum(Data)
phi <- extraDistr::rdirichlet(k, alpha0*1000)

#混合率の初期値
r <- rep(1/k, k)

#観測データの対数尤度の初期化
L <- LLobz(Data, phi, r, const, hh, k)
LL1 <- L$LLob
z <- L$z

#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.5
iter <- 0 

##EMアルゴリズムで対数尤度を最大化
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #Eステップの計算
  z <- L$z   #潜在変数zの出力
  
  #Mステップの計算と最適化
  #phiの推定
  phi <- matrix(0, nrow=k, ncol=v)
  for(j in 1:k){
    #完全データの対数尤度からphiの推定量を計算
    phi[j, ] <- colSums(matrix(z[, j], nrow=hh, ncol=v) * Data) / sum(z[, j] * s)   #重み付き多項分布の最尤推定
  }
  
  #混合率を推定
  r <- apply(z, 2, sum) / hh
  
  #観測データの対数尤度を計算
  phi[phi==0] <- min(phi[phi > 0])
  L <- LLobz(Data, phi, r, const, hh, k)
  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####混合多項分布モデルの推定結果####
##EMアルゴリズムの推定結果と真値の比較
round(cbind(t(phi), t(phit)), 3)   #phiの推定値とphiの真値
round(cbind(Z %*% 1:k, apply(z, 1, which.max), z), 3)   #潜在変数zとセグメント割当
round(r, 3)   #混合率

##適合度の確認
L$LLob   #観測データの対数尤度
-2*(L$LLob) + 2*length(phi)   #AIC



