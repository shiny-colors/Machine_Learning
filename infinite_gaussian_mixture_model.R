#####無限混合ガウス分布モデル#####
library(MASS)
library(mclust)
library(reshape2)
library(gtools)
library(bayesm)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

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
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####データの発生####
##データの設定
hh <- 5000
seg <- 5
k <- 4

#セグメント割当の設定
seg_id <- rep(1:seg, rep(hh/seg, seg))

##パラメータの設定
#平均構造の発生
mu0 <- matrix(runif(k*seg, -5, 5), nrow=seg, ncol=k)

#分散共分散行列の発生
Cov0 <- list()
for(j in 1:seg){
  Cor <- corrM(k, -0.55, 0.7, 0.1, 0.3)
  Cov0[[j]] <- covmatrix(k, Cor, 0.5, 2)$covariance
}

##多変量正規分布からデータを発生
Data <- matrix(0, nrow=hh, ncol=k)
for(j in 1:seg){
  Data[seg_id==j, ] <- mvrnorm(length(seg_id[seg_id==j]), mu0[j, ], Cov0[[j]])
}

#二変量ごとの結果を可視化
plot(Data[, 1:2], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ1の値", ylab="データ2の値", main="混合二変量正規分布のプロット")
plot(Data[, c(1, 3)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ1の値", ylab="データ3の値", main="混合二変量正規分布のプロット")
plot(Data[, c(1, 4)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ1の値", ylab="データ4の値", main="混合二変量正規分布のプロット")
plot(Data[, 2:3], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ2の値", ylab="データ3の値", main="混合二変量正規分布のプロット")
plot(Data[, c(2, 4)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ2の値", ylab="データ4の値", main="混合二変量正規分布のプロット")
plot(Data[, 3:4], col=seg_id[seg_id %in% 1:5], pch=20, xlab="データ3の値", ylab="データ4の値", main="混合二変量正規分布のプロット")


####マルコフ連鎖モンテカルロ法で無限混合ガウス分布モデルを推定####
##アルゴリズムの設定
R <- 10000
keep <- 2
sbeta <- 1.5
iter <- 0

##多変量正規分布の尤度関数
dmv <- function(x, mean.vec, S, S_det, S_inv){
  LLo <- (2*pi)^(-nrow(S)/2) * S_det^(-1/2) *
    exp(-1/2 * (x - mean.vec) %*% S_inv %*% (x - mean.vec))
  return(LLo)
}

##事前分布の設定
#多変量正規分布の事前分布
mean0 <- rep(0, k)
sigma0 <- diag(100, k)
sigma0_inv <- solve(sigma0)
mean_z <- colMeans(Data)   #提案分布の平均
sigma_z <- diag(diag(var(Data))) * 7.5   #提案分布の分散

#逆ウィシャート分布の事前分布
nu <- k+1
V <- nu * diag(k)

#ディクレリ分布の事前分布
alpha <- 1

##初期値の設定
seg0 <- 2   #初期セグメントは2つ
res <- Mclust(Data, 2)   #混合正規分布を推定

#潜在変数zの初期値
z_vec <- apply(res$z, 1, which.max)
z <- matrix(0, nrow=hh, ncol=seg0)
for(i in 1:hh) {z[i, z_vec[i]] <- 1}

#回帰パラメータの初期値
oldmean <- t(res$parameters$mean)
oldcov <- res$parameters$variance$sigma

##パラメータの格納用配列
max_seg <- 15
Z <- matrix(0, nrow=R/keep, ncol=max_seg)
Mu <- array(0, dim=c(max_seg, k, R/keep))
Cov <- list()
storage.mode(Z) <- "integer"


####MCMCでパラメータをサンプリング
for(rp in 1:R){
  
  ##潜在変数zをサンプリング
  z_len <- length(unique(z_vec))
  LLi <- matrix(0, nrow=hh, ncol=z_len+1)
  
  #サンプリング済みの潜在変数の尤度
  for(j in 1:z_len){
    LLi[, j] <- dmvnorm(Data, oldmean[j, ], oldcov[, , j])
  }
  #新しい潜在変数の尤度
  LLi[, z_len+1] <- dmvnorm(Data, mean_z, sigma_z)
  
  #CRPを計算
  gamma0 <- cbind(matrix(colSums(z), nrow=hh, ncol=z_len, byrow=T) - z, alpha)
  gamma1 <- LLi * gamma0/(hh-1-alpha)
  
  #多項分布より潜在変数zをサンプリング
  z <- t(apply(gamma1/rowSums(gamma1), 1, function(x) rmultinom(1, 1, x)))
  z <- z[, colSums(z) > 0]
  z_vec <- z %*% 1:ncol(z)
  
  ##多変量正規分布のパラメータをサンプリング
  #生成されなかったzは消去しておく
  z_cnt <- length(colSums(z) > 0)
  if(z_cnt > nrow(oldmean)){
    oldmean <- rbind(oldmean, 0)
  }
  oldmean <- oldmean[1:z_cnt, ]
  oldcov <- array(0, dim=c(k, k, z_cnt))
  
  #多変量正規分布の平均および分散共分散行列をギブスサンプリング
  for(j in 1:z_cnt){
    
    #インデックスを作成
    index <- subset(1:nrow(z), z[, j]==1)
    
    #逆ウィシャート分布から分散共分散行列をサンプリング
    Vn <- nu + length(index)
    er <- Data[index, ] - matrix(oldmean[j, ], nrow=length(index), ncol=ncol(Data), byrow=T)
    R_par <- solve(V) + t(er) %*% er
    
    oldcov[, , j] <- rwishart(Vn, solve(R_par))$IW   #逆ウィシャート分布から分散共分散行列をサンプリング
    
    #平均をサンプリング
    if(length(index) > 1){
      mean_mu <- length(index)/(1+length(index))*colMeans(Data[index, ])   #平均パラメータの平均
    } else {
      mean_mu <- length(index)/(1+length(index))*Data[index, ] 
    }
    mean_cov <- oldcov[, , j] / (1+length(index))   #平均パラメータの分散共分散行列
    oldmean[j, ] <- mvrnorm(1, mean_mu, mean_cov)
  }
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    Mu[1:nrow(oldmean), , mkeep] <- oldmean
    Cov[[mkeep]] <- oldcov
    if(rp >= R/2){Z[, 1:ncol(z)] <- Z[, 1:ncol(z)] + z}   #繰り返し数が最大反復数の半分を超えたらパラメータを格納
    
    print(rp)
    print(colSums(z))
  }
}

####サンプリング結果の要約と可視化####
#バーンイン期間
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##サンプリング結果をプロット
matplot(t(Mu[1, , burnin2:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(Mu[2, , burnin2:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(Mu[3, , burnin2:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(Mu[4, , burnin2:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(Mu[5, , burnin2:(R/keep)]), type="l", ylab="パラメータ")

##サンプリング結果の事後平均
mcmc_seg <- sum(colSums(Z) > 10000)   #推定されたセグメント数

#潜在変数zの推定量
round(Z_mu <- (Z/rowSums(Z))[, colSums(Z) > 0], 3)   #潜在変数の割当確率
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #混合率

#平均の推定量
mean_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  mean_mu[i, ] <- colMeans(t(Mu[i, , burnin1:(R/keep)]))
}
round(rbind(mean_mu, mu0), 3)   #真のパラメータと比較

#分散共分散行列の推定量
cov_mu0 <- Cov[[burnin1]]
for(i in (burnin1+1):(R/keep)){
  cov_mu0 <- cov_mu0 + Cov[[i]][, , 1:mcmc_seg]
}
round(Cov_mu <- cov_mu0/length(burnin1:(R/keep)), 3)   #分散共分散行列の事後平均
lapply(Cov0, function(x) round(x, 3))   #真の分散共分散行列

