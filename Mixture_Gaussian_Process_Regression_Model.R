######混合ガウス過程回帰モデル#####
options(warn=2)
library(MASS)
library(kernlab)
library(GPfit)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5698)

####データの発生####
#データの設定
seg <- 8   #混合数
d <- 5000   #サンプル数
k <- 100   #入力変数数
w <- rpois(d, rgamma(d, 15.0, 0.5))
w[w < 3] <- ceiling(runif(sum(w < 3), 3, 10))

#セグメントの生成
alpha01 <- as.numeric(extraDistr::rdirichlet(1, rep(20.0, seg)))
Z <- rmnom(d, 1, alpha01)
z <- as.numeric(Z %*% 1:seg)

index_z <- list()
for(j in 1:seg){
  index_z[[j]] <- which(z==j)
} 


##データの生成
#多項分布から入力変数を生成
for(j in 1:1000){
  alpha0 <- rep(0.15, k)
  theta <- thetat <- extraDistr::rdirichlet(seg, alpha0)   #多項分布のパラメータ
  Data <- rmnom(d, w, theta[z, ])
  if(min(colSums(Data)) >= 10) break
}
sparse_data <- as(Data, "CsparseMatrix")


#カーネル関数の生成
kern_data <- list()
n0 <- rep(0, seg)
for(j in 1:seg){
  index <- index_z[[j]]
  kern_data[[j]] <- Data[index, ] %*% t(Data[index, ])
  n0[j] <- nrow(kern_data[[j]])
}

#ガウス過程からセグメントごとに応答変数を生成
y <- rep(0, d)
sigma <- rep(0, seg)
for(j in 1:seg){
  K <- kern_data[[j]]
  n <- n0[j]
  sigma[j] <- runif(1, 0.5, 5.0)
  y[index_z[[j]]] <- mvrnorm(1, rep(0, n0[j]), K + diag(sigma[j], n))   #ガウス過程から応答変数を生成
}


####MCMC-EMアルゴリズムで混合ガウス過程回帰モデルを推定####
##アルゴリズムの設定
R <- 5000   #サンプリング回数
keep <- 2   #2回に1回の割合でサンプリング結果を格納
disp <- 10
iter <- 0
burnin <- 1000/keep

#事前分布の設定
alpha01 <- rep(10, seg)
alpha02 <- rep(1, k)
r <- matrix(1/seg, nrow=d, ncol=seg, byrow=T)

#初期値の設定
theta <- extraDistr::rdirichlet(seg, rep(5, k))
tau1 <- rep(1, seg)
tau2 <- rep(0.5, seg)

#パラメータの格納用配列
THETA <- array(0, dim=c(seg, k, R/keep))
SEG <- matrix(0, nrow=d, ncol=seg)
TAU1 <- matrix(0, nrow=R/keep, ncol=seg)
TAU2 <- matrix(0, nrow=R/keep, ncol=seg)

#データの設定
const <- lfactorial(w) - rowSums(lfactorial(Data))


####パラメータをサンプリング####
for(rp in 1:R){
  
  ##混合多項分布モデルより入力変数の割当をサンプリング
  #多項分布の対数尤度関数
  LLi <- as.matrix(const + sparse_data %*% t(log(theta)))

  #潜在変数の割当Zをサンプリング
  LL_max <- rowMaxs(LLi)
  LH <- exp(LLi - LL_max)   #尤度に変換
  z_rate <- r * LH / rowSums(r * LH)   #潜在変数の割当確率
  Zi <- rmnom(d, 1, z_rate)   #潜在変数をサンプリング
  z_vec <- as.numeric(Zi %*% 1:seg)
  
  ##混合多項分布のパラメータを更新
  #混合率の更新
  rsum <- colSums(Zi) + alpha01
  r <- matrix(extraDistr::rdirichlet(1, rsum), nrow=d, ncol=seg, byrow=T)

  #多項分布のパラメータを更新
  wsum0 <- t(Data) %*% Zi
  wsum <- t(wsum0 + alpha02)
  theta <- extraDistr::rdirichlet(seg, wsum)   #ディリクレ分布からthetaをサンプリング
  
  
  ##セグメント割当に基づきガウス過程回帰モデルをEMアルゴリズムで推定
  index_zi <- list()
  n <- rep(0, seg)
  LLm <- rep(0, j)
  beta_list <- list()
  
  for(j in 1:seg){
    #潜在変数zのインデックスを作成
    index_zi[[j]] <- which(Zi[, j]==1)
    index <- index_zi[[j]]
  
    #データとカーネル関数を設定
    y_vec <- y[index]
    data <- Data[index, ]
    K <- data %*% t(data)
    KK <- K %*% K
    n[j] <- length(index_zi[[j]])
    n0 <- n[j]
    tau_s1 <- tau1[j]
    tau_s2 <- tau2[j]
  
    #Eステップで回帰係数を推定
    beta <- solve(KK + diag(tau_s1/tau_s2, n0)) %*% t(K) %*% y_vec
    delta <- diag(tau_s1, n0) + tau_s2*KK
    beta_list[[j]] <- beta
    
    #Mステップでハイパーパラメータを推定
    tau1_inv <- (sum(abs(beta^2)) + sum(diag(solve(delta)))) / n0
    tau1[j] <- 1 / tau1_inv
    tau2_inv <- (sum((y_vec - K %*% beta)^2) + sum(diag(KK %*% solve(delta)))) / n0
    tau2[j] <- 1 / tau2_inv
    
    #周辺尤度の更新
    diag_tau2 <- diag(tau2_inv, n0)
    tau1_KK <- tau1_inv * KK
    Lm <- -1/2*abs(diag_tau2 + tau1_KK) - 1/2*as.numeric((t(y_vec) %*% solve(diag_tau2 + tau1_KK) %*% y_vec))
    LLm[j] <- sum(Lm)
  }
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    TAU1[mkeep, ] <- tau1
    TAU2[mkeep, ] <- tau2
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(rp%%keep==0 & rp >= burnin){
      SEG <- SEG + Zi
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(LLm)
      print(round(rbind(r[1, ], colMeans(Z)), 3))
      print(round(cbind(theta[, 1:10], thetat[, 1:10]), 3))
    }
  }
}

####推定結果の確認と要約####
burnin <- 1000/keep
RS <- R/keep

##サンプリング結果の可視化
#多項分布のパラメータの可視化
matplot(t(THETA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[7, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

#正則化パラメータの可視化
matplot(TAU1/TAU2, type="l", xlab="サンプリング回数", ylab="パラメータ")


##推定結果の要約推定量
#正則化パラメータの要約統計量
tau_mu <- colMeans(TAU1[burnin:RS, ])/colMeans(TAU2[burnin:RS, ])

##潜在変数の割当ごとに回帰係数を推定
#潜在変数の割当
Zi <- SEG / rowSums(SEG)

#重み付きデータを設定
data_list <- list()
K_list <- list()
y_list <- list()
beta_list <- list()

for(j in 1:seg){
  data0 <- Zi[, j] * Data; y_vec0 <- Zi[, j] * y   
  data_list[[j]] <- data <- data0[rowSums(abs(data0)) > 0, ]
  y_list[[j]] <- y_vec <- y_vec0[abs(y_vec0) > 0]
  n <- nrow(data)
  
  #回帰係数を推定
  K_list[[j]] <- K <- data %*% t(data)   #カーネル関数
  KK <- K %*% K
  beta_list[[j]] <- solve(KK + diag(tau_mu, n)) %*% t(K) %*% y_vec
}

#予測結果を推定
matplot(cbind(K_list[[1]] %*% beta_list[[1]], y_list[[1]]), type="l", xlab="サンプル番号", ylab="パラメータ")
matplot(cbind(K_list[[3]] %*% beta_list[[3]], y_list[[3]]), type="l", xlab="サンプル番号", ylab="パラメータ")
matplot(cbind(K_list[[5]] %*% beta_list[[5]], y_list[[5]]), type="l", xlab="サンプル番号", ylab="パラメータ")
matplot(cbind(K_list[[7]] %*% beta_list[[7]], y_list[[7]]), type="l", xlab="サンプル番号", ylab="パラメータ")

