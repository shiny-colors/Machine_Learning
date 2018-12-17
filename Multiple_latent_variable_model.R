#####Multiple lattent variable model#####
library(MASS)
library(matrixStats)
library(Matrix)
library(data.table)
library(bayesm)
library(extraDistr)
library(condMVNorm)
library(gtools)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(78594)

####任意の分散共分散行列を作成させる関数####
##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  X.Sigma <- eigen(Sigma)
  Lambda <- diag(X.Sigma$values)
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda)
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  diag(Sigma) <- 1
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
k <- 15   #クラスタ数
hh <- 5000   #ユーザー数
pt <- rtpois(hh, rgamma(hh, 10.0, 0.55), a=0, b=Inf)   #ユーザーあたりのレコード数
hhpt <- sum(pt)   #総レコード数

##IDとインデックスの設定
#IDの設定
no <- 1:hhpt
u_id <- rep(1:hh, pt)
t_id <- as.numeric(unlist(tapply(1:hhpt, u_id, rank)))

#インデックスの設定
user_list <- list()
for(i in 1:hh){
  user_list[[i]] <- which(u_id==i)
}
u_dt <- sparseMatrix(1:hhpt, u_id, x=rep(1, hhpt), dims=c(hhpt, hh))


##応答変数が妥当な数値になるまで繰り返す
rp <- 0
repeat {
  rp <- rp + 1

  ##パラメータの設定
  #多変量正規分布のパラメータを設定
  mut <- mu <- rnorm(k, -0.2, 0.5)
  Covt <- Cov <- corrM(k, -0.7, 0.8, 0.05, 0.2)
  
  #回帰ベクトルを生成
  betat <- beta <- rnorm(k, -0.4, 1.5)
  
  ##応答変数を生成
  #多変量正規分布からクラスタを生成
  U <- mvrnorm(hh, mut, Cov)   
  Z <- matrix(as.numeric(U > 0), nrow=hh, ncol=k)
  mean(Z)
  
  #ロジットと選択確率を設定
  logit <- as.numeric(Z[u_id, ] %*% beta)
  Prob <- exp(logit) / (1 + exp(logit))
  
  #ベルヌーイ分布から応答変数を生成
  y <- rbinom(hhpt, 1, Prob)
  print(mean(y))
  
  if(mean(y) > 0.2 & mean(y) < 0.4){
    break
  }
}

####マルコフ連鎖モンテカルロ法でMultiple lattent variable modelを推定####
##切断正規分布の乱数を発生させる関数
rtnorm <- function(mu, sigma, a, b){
  FA <- pnorm(a, mu, sigma)
  FB <- pnorm(b, mu, sigma)
  return(qnorm(runif(length(mu))*(FB-FA)+FA, mu, sigma))
}

##多変量正規分布の条件付き期待値と条件付き分散を計算する関数
cdMVN <- function(mean, Cov, dependent, U){
  
  #分散共分散行列のブロック行列を定義
  Cov11 <- Cov[dependent, dependent]
  Cov12 <- Cov[dependent, -dependent, drop=FALSE]
  Cov21 <- Cov[-dependent, dependent, drop=FALSE]
  Cov22 <- Cov[-dependent, -dependent]
  
  
  #条件付き分散と条件付き平均を計算
  CDinv <- Cov12 %*% solve(Cov22)
  CDmu <- mean[, dependent] + t(CDinv %*% t(U[, -dependent] - mean[, -dependent]))   #条件付き平均を計算
  CDvar <- Cov11 - Cov12 %*% solve(Cov22) %*% Cov21   #条件付き分散を計算
  val <- list(CDmu=CDmu, CDvar=CDvar)
  return(val)
}

##多変量正規分布の密度関数
mvdnorm <- function(u, mu, Cov, s){
  er <- u - mu   #誤差
  Lho <- 1 / (sqrt(2*pi)^s*sqrt(det(Cov))) * exp(-1/2 * as.numeric((er %*% solve(Cov) * er) %*% rep(1, s)))
  return(Lho)
}

##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10


Z





