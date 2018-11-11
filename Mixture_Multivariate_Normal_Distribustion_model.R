#####高速混合多変量正規分布#####
library(MASS)
library(matrixStats)
library(FAdist)
library(mclust)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
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
k <- 12   #混合数
N <- 100000   #サンプル数
v <- 7   #次元数
vec <- rep(1, k)

##パラメータの設定
#平均ベクトルを生成
mu <- matrix(0, nrow=k, ncol=v)
for(j in 1:k){
  lower <- runif(v, 30, 70)   #変数の下限値
  upper <- lower + runif(v, 30, 60)   #変数の上限値
  mu[j, ] <- runif(v, lower, upper)   #変数の平均値を発生
}
mut <- mu

#分散共分散行列を生成
Cov <- array(0, dim=c(v, v, k))
lower <- 4; upper <- 25
for(j in 1:k){
  Cor <- corrM(v, -0.6, 0.9, 0.1, 0.5)
  Cov[, , j] <- covmatrix(v, Cor, lower, upper)$covariance
}
Covt <- Cov

##応答変数を生成
#セグメント割当を生成
prob <- extraDistr::rdirichlet(N, rep(25.0, k))
Z <- rmnom(N, 1, prob)
z_vec <- as.numeric(Z %*% 1:k)
r <- colMeans(Z)

#多変量正規分布からデータを生成
y <- matrix(0, nrow=N, ncol=v)
for(i in 1:N){
  y[i, ] <- mvrnorm(1, mu[z_vec[i], ], Cov[, , z_vec[i]])
}


####EMアルゴリズムで混合多変量正規分布を推定####
##多変量正規分布の密度関数
mvdnorm <- function(u, mu, Cov, N, s){
  er <- y - matrix(mu, nrow=N, ncol=v, byrow=T)
  Lho <- 1 / (sqrt(2*pi)^s*sqrt(abs(det(Cov)))) * exp(-1/2 * as.numeric((er %*% ginv(Cov) * er) %*% rep(1, s)))
  return(Lho)
}

##観測データの対数尤度と潜在変数zの定義
LLobz <- function(y, mu, Cov, r, v, k, vec){
  #混合多変量正規分布のセグメントごとの尤度
  LLind <- matrix(0, nrow=N, ncol=k)
  for(j in 1:k){
    LLind[, j] <- mvdnorm(y, mu[j, ], Cov[, , j], N, v)   #多変量正規分布の尤度
  }

  #対数尤度と潜在変数zを定義
  LLho <- matrix(r, nrow=N, ncol=k, byrow=T) * LLind; LLho_vec <- as.numeric(LLho %*% vec)
  z <- LLho / LLho_vec   #潜在変数zの割当確率
  LL <- sum(log(LLho_vec))   #観測データの対数尤度の和
  value <- list(LL=LL, z=z)
  return(value)
}

##EMアルゴリズムの設定
dl <- 100   #EMステップでの対数尤度の差の初期化
tol <- 0.001 
iter <- 1

##パラメータの初期値
Zi <- rmnom(N, 1, rep(20, k))
mu <- matrix(0, nrow=k, ncol=v)
Cov <- array(0, dim=c(v, v, k))
r <- rep(1/k, k)
for(j in 1:k){
  data <- y[Zi[, j]==1, ]
  mu[j, ] <- colMeans(data) + runif(v, -15, 15)
  Cov[, , j] <- var(data)
}

#対数尤度の初期化
Lobz <- LLobz(y, mu, Cov, r, v, k, vec)
LL1 <- Lobz$LL


####EMアルゴリズムによるパラメータの更新####
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  ##Mステップでパラメータを最尤推定
  z <- Lobz$z   #潜在変数zの出力
  
  for(j in 1:k){
    #平均ベクトルを最尤推定
    weighted_z <- sum(z[, j])
    weighted_data <- z[, j] * y   #重み付きの応答変数
    mu[j, ] <- colSums(weighted_data) / weighted_z
    
    #分散共分散行列を最尤推定
    er <- weighted_data - (z[, j] * matrix(mu[j, ], nrow=N, ncol=v, byrow=T))   #重み付き誤差
    Cov[, , j] <- t(er) %*% er / weighted_z
  }
  #混合率を更新
  r <- colSums(z) / N
  
  ##Eステップで観測データの対数尤度と潜在変数zの定義
  Lobz <- LLobz(y, mu, Cov, r, v, k, vec)   #観測データの対数尤度
  LL <- Lobz$LL
  
  ##アルゴリズムの更新
  ite <- iter + 1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}
res <- Mclust(y, k)

