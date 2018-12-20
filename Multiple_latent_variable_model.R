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
u_dt <- sparseMatrix(u_id, 1:hhpt, x=rep(1, hhpt), dims=c(hh, hhpt))


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
  UT <- U <- mvrnorm(hh, mut, Cov)   
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

##テストデータを作成
#ロジットと確率を設定
test_id <- rep(1:hh, rtpois(hh, 10.0, a=0, b=Inf))
f <- length(test_id)
logit <- as.numeric(Z[test_id, ] %*% beta)
Prob <- exp(logit) / (1 + exp(logit))

#ベルヌーイ分布から応答変数を生成
y_test <- rbinom(f, 1, Prob)

 

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

##リープフロッグ法を解く関数
leapfrog <- function(r, z, D, e, L) {
  leapfrog.step <- function(r, z, e){
    r2 <- r  - e * D(z, Zi_dt, y) / 2
    z2 <- z + e * r2
    r2 <- r2 - e * D(z2, Zi_dt, y) / 2
    list(r=r2, z=z2) # 1回の移動後の運動量と座標
  }
  leapfrog.result <- list(r=r, z=z)
  for(i in 1:L) {
    leapfrog.result <- leapfrog.step(leapfrog.result$r, leapfrog.result$z, e)
  }
  leapfrog.result
}

##ロジスティック回帰モデルの対数尤度関数と対数微分関数
#ロジスティック回帰モデルの対数尤度関数を定義
loglike <- function(beta, Data, y){
  
  #ロジットと確率を定義
  mu <- exp(as.numeric(Data %*% beta))
  prob <- mu / (1 + mu)
  
  #対数尤度の和
  LL <- sum(y*log(prob) + (1-y)*log(1-prob))  
  return(LL)
}

#ロジスティック回帰モデルの対数尤度の微分関数
dloglike <- function(beta, Data, y){
  
  #ロジットと確率を定義
  mu <- exp(as.numeric(Data %*% beta))
  prob <- mu / (1 + mu)
  
  #勾配ベクトルを計算
  dlogit <- y*Data - Data*prob
  LLd <- -colSums(dlogit)
  return(LLd)
}

##アルゴリズムの設定
R <- 2500
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10
e <- 0.0025
L <- 3

##事前分布の設定
#回帰ベクトルの事前分布
alpha1 <- rep(0, k)

#多変量正規分布の事前分布
alpha2 <- rep(0, k)   #平均ベクトルの事前分布
inv_tau <- solve(diag(100, k))
nu <- 1   #逆ウィシャート分布の自由度
V <- 2.0*k * diag(k)   #逆ウィシャート分布のパラメータ

##パラメータの真値
#モデルパラメータの真値
beta <- betat
mu <- mut
mu_dt <- matrix(mu, nrow=hh, ncol=k, byrow=T)
Cov <- Covt

#潜在変数の真値
U <- UT
Zi <- Z

##パラメータの初期値
#モデルパラメータの初期値
beta <- runif(k, -0.5, 0.5)
mu <- rep(-0.25, k)
mu_dt <- matrix(mu, nrow=hh, ncol=k, byrow=T)
Cov <- diag(1, k)

#潜在変数の初期値
U <- mvrnorm(hh, mu, Cov)
Zi <- matrix(as.numeric(U > 0), nrow=hh, ncol=k)

##サンプリング結果の格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=k)
MU <- matrix(0, nrow=R/keep, ncol=k)
COV <- array(0, dim=c(k, k, R/keep))
SEG <- matrix(0, nrow=hh, ncol=k)

##切断正規分布の切断領域を定義
a <- ifelse(Zi==0, -100, 0)
b <- ifelse(Zi==1, 100, 0)

##対数尤度の基準値
#1パラメータモデルの対数尤度
LLst <- sum(y_test*log(mean(y_test))) + sum((1-y_test)*log(1-mean(y_test)))

#ベストなパラメータでの対数尤度
logit <- as.numeric(Z[test_id, ] %*% betat)
prob <- exp(logit) / (1 + exp(logit))
LLbest <- sum(y_test*log(prob) + (1-y_test)*log(1-prob))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##潜在効用からクラスタ割当をサンプリング
  for(j in 1:k){
    
    #多変量正規分布の条件付き分布からクラスタ割当の事前確率を設定
    MVR <- cdMVN(matrix(mu_dt, nrow=hh, ncol=k, byrow=T), Cov, j, U)
    MVR_U <- as.numeric(MVR$CDmu); MVR_S <- sqrt(MVR$CDvar)
    prob_z <- pnorm(MVR_U, 0, MVR_S)   #分布関数から事前確率を設定
    
    #クラス割当ごとのロジットと応答確率を設定
    logit1 <- as.numeric(Zi[u_id, -j] %*% beta[-j])
    logit2 <- logit1 + beta[j] 
    prob1 <- exp(logit1) / (1 + exp(logit1))
    prob2 <- exp(logit2) / (1 + exp(logit2))
    
    #ユーザーごとの尤度の積
    Li1 <- exp(as.numeric(u_dt %*% (y*log(prob1) + (1-y)*log(1-prob1))))
    Li2 <- exp(as.numeric(u_dt %*% (y*log(prob2) + (1-y)*log(1-prob2))))
    
    #ベルヌーイ分布からクラスタ割当を生成
    prob_class <- (Li2*prob_z) / (Li1*(1-prob_z) + Li2*prob_z)   #ベイズの定理からクラスタ割当確率を計算
    Zi[, j] <- rbinom(hh, 1, prob_class)
    
    #切断正規分布から潜在効用をサンプリング
    a[, j] <- ifelse(Zi[, j]==0, -100, 0)
    b[, j] <- ifelse(Zi[, j]==1, 100, 0)
    U[, j] <- rtnorm(MVR_U, MVR_S, a[, j], b[, j])
    U[is.infinite(U[, j])==TRUE, j] <- 0
  }
  
  ##HMCで回帰ベクトルをサンプリング
  #パラメータを設定
  Zi_dt <- Zi[u_id, ]
  rold <- rnorm(k)
  betad <- beta
 
  res <- leapfrog(rold, betad, dloglike, e, L)   #リープフロッグ法による1ステップ移動
  rnew <- res$r
  betan <- res$z

  #移動前と移動後のハミルトニアンを計算
  Hnew <- -loglike(betan, Zi_dt, y) + sum(rnew^2)/2
  Hold <- -loglike(betad, Zi_dt, y) + sum(rold^2)/2
  
  #HMC法によりパラメータの採択を決定
  gamma <- min(1, exp(Hold - Hnew))
  if(gamma=="NaN") gamma <- -1
  
  #一様乱数を発生
  u <- runif(1)
  
  #u < alphaなら新しいbetaを採択
  if(u < gamma){
    beta <- betan
    #そうでないならbetaを更新しない
  } else {
    beta <- betad
  }
  
  ##多変量正規分布のパラメータをサンプリング
  #平均ベクトルをサンプリング
  mu_vec <- colMeans(U)
  inv_Sigma <- solve(inv_tau + diag(hh, k))
  mu_par <- as.numeric(inv_Sigma %*% (diag(hh, k) %*% mu_vec))
  mu <- mvrnorm(1, mu_par, inv_Sigma)   #多変量正規分布から平均ベクトルをサンプリング
  
  #相関行列をサンプリング
  #逆ウィシャート分布のパラメータ
  R_error <- U - matrix(mu_vec, nrow=hh, ncol=k, byrow=T) 
  IW_R <- V + t(R_error) %*% R_error   #逆ウィシャート分布のパラメータ
  Sn <-  hh + nu   #逆ウィシャート分布の自由度
  Cov_hat <- rwishart(Sn, solve(IW_R))$IW   #逆ウィシャート分布からパラメータをサンプリング
  Cov <- cov2cor(Cov_hat)
  
  
  ##パラメータの格納とサンプリング結果の表示
  #パラメータを格納
  if(rp%%keep==0){
    #モデルのパラメータを格納
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    MU[mkeep, ] <- mu
    COV[, , mkeep] <- Cov
    
    #バーンイン期間を超えたらトピックを格納
    if(rp >= burnin){
      SEG <- SEG + Zi
    }
  }
  
  #対数尤度の計算とサンプリング結果の表示
  if(rp==1 | rp%%disp==0){
    #対数尤度を計算
    logit <- as.numeric(Zi[test_id, ] %*% beta)
    prob <- exp(logit) / (1 + exp(logit))
    LL <- sum(y_test*log(prob) + (1-y_test)*log(1-prob))
    
    #サンプリング結果を表示
    print(rp)
    print(gamma)
    print(round(Cov, 3))
    print(c(LL, LLst, LLbest))
  }
}

