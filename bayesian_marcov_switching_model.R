#####ベイジアン一般状態空間マルコフ切換えモデル#####
library(MASS)
library(MSwM) 
library(matrixStats)
library(FAdist)
library(Rfast)
library(bayesm)
library(extraDistr)
library(actuar)
library(gtools)
library(caret)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(18478)

####データの発生####
n <- 3000   #観測期間数
s <- 4   #切換え数
k1 <- 2   #観測モデルの応答変数のパラメータ数
k2 <- 2   #観測モデルの説明変数のパラメータ数

####パラメータの設定####
##マルコフ切換え行列のパラメータを設定
#初期確率の設定
p0 <- c(0.5, 0.2, 0.2, 0.1)

#推移行列の設定
pr1 <- c(0.7, 0.3, 0, 0)
pr2 <- c(0, 0, 0.6, 0.4)
pr3 <- c(0.75, 0.25, 0, 0)
pr4 <- c(0, 0, 0.8, 0.2)
Pr0 <- rbind(pr1, pr2, pr3, pr4)

##観測モデルのパラメータを設定
#標準化価格のパラメータ
tau01 <- c(8.0, 7.0, 9.0, 7.5)
tau02 <- c(0.5, 3.5, 1.0, 2.25)
tau0 <- cbind(tau01, tau02)
colnames(tau0) <- c("", "")


#売上点数のパラメータ
beta00 <- c(2.1, 3.3, 1.5, 2.7)
beta01 <- c(1.0, 1.8, 0.8, 1.3)
beta02 <- c(0.5, 0.7, 0.3, 0.6)
beta0 <- rbind(beta00, beta01, beta02)
sigma0 <- 0.5


####応答変数の発生####
#データの保存用配列
Z0 <- matrix(0, nrow=n, ncol=s)
z0 <- rep(0, n)
y <- matrix(0, nrow=n, ncol=k1)
X <- cbind(1, matrix(0, nrow=n, ncol=k2))
X[, ncol(X)] <- rbinom(n, 1, 0.4)

#初期値の設定
Z0[1, ] <- extraDistr::rmnom(1, 1, p0)
seg <- z0[1] <- Z0[1, ] %*% 1:s
y[1, 2] <- rbeta(1, tau01[seg], tau02[seg])
X[1, 2] <- y[1, 2]
y[1, 1] <- X[1, ] %*% beta0[, seg] + rnorm(1, 0, sigma0)

##時間ごとに逐次的に応答変数を発生させる
for(i in 2:n){
  print(i)
  #マルコフ推移行列から潜在変数の割当を生成
  obs_seg <- z0[i-1]
  Z0[i, ] <- extraDistr::rmnom(1, 1, Pr0[obs_seg, ])
  seg <- z0[i] <- Z0[i, ] %*% 1:s
  
  #潜在変数から観測データを発生
  y[i, 2] <- rbeta(1, tau01[seg], tau02[seg])
  X[i, 2] <- y[i, 2]
  y[i, 1] <- X[i, ] %*% beta0[, seg] + rnorm(1, 0, sigma0)
}
table(z0)


####マルコフ連鎖モンテカルロ法でマルコフスイッチングモデルを推定####
##アルゴリズムの設定
R <- 20000  
keep <- 4
iter <- 0
sbeta <- 1.5

##事前分布の設定
#線形回帰モデルの事前分布
alpha0 <- rep(0, k2+1)
cov0 <- diag(0.01, k2+1)
s0 <- 0.01
v0 <- 0.01

#ベータ分布モデルの事前分布
lambda0 <- rep(0, 2)
omega0 <- diag(0.01, 2)

#マルコフ推移行列の事前分布
phi01 <- rep(1, 2)
phi02 <- rep(1, s)

##初期値の設定
#線形回帰モデルの初期値
oldbeta <- matrix(solve(t(X) %*% X) %*% t(X) %*% y[, 1], nrow=k1+1, ncol=s)
oldbeta[1, ] <- c(1.5, 2.5, 1, 2)
oldsigma <- sd(rep(y[, 1], s) - as.numeric(X %*% oldbeta))

#ベータ分布の初期値
oldtau <- matrix(beta.mle(y[, 2])$param, nrow=s, ncol=2, byrow=T)
oldtau[, 2] <- c(0.5, 2.0, 0.75, 1.5)
rw <- rbind(c(0.005, 0.0025), matrix(c(0.01, 0.01), nrow=s-1, ncol=2, byrow=T))


#推移行列の初期値
Pr1 <- c(0.7, 0.3, 0, 0)
Pr2 <- c(0, 0, 0.7, 0.3)
Pr3 <- c(0.7, 0.3, 0, 0)
Pr4 <- c(0, 0, 0.7, 0.3)
Pr <- rbind(Pr1, Pr2, Pr3, Pr4)
r0 <- r <- rep(0.25, s)

##インデックスの設定
index_s <- list()
for(i in 1:s){
  index_s[[i]] <- which(Pr[i, ] > 0)
}

##パラメータの格納用配列
BETA <- matrix(0, nrow=R/keep, ncol=(k1+1)*s)
SIGMA <- rep(0, R/keep)
LAMBDA <- matrix(0, nrow=R/keep, ncol=k2*s)
Zi <- matrix(0, nrow=R/keep, ncol=n)
PR <- array(0, dim=c(s, s, R/keep))
storage.mode(Zi) <- "integer"


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##時間ごとに動的な潜在変数zを生成(システムモデルのパラメータをサンプリング)
  #データの格納用配列
  Z <- matrix(0, nrow=n, ncol=s)
  z <- rep(0, n)
  z_rate <- matrix(0, nrow=n, ncol=s)
  LLi1 <- matrix(0, nrow=n, ncol=s)
  LLi2 <- matrix(0, nrow=n, ncol=s)
  
  #観測モデルの尤度を計算
  mu <- X %*% oldbeta   #線形回帰モデルの平均構造
  for(j in 1:s){
    LLi1[, j] <- dnorm(y[, 1], mu[, j], oldsigma)
    LLi2[, j] <- dbeta(y[, 2], oldtau[j, 1], oldtau[j, 2])
  }
  LLi <- LLi1 * LLi2
  
  #多項分布より1期目の潜在変数zを生成
  z_rate[1, ] <- r*LLi[1, ] / sum(r*LLi[1, ])
  Z[1, ] <- extraDistr::rmnom(1, 1, z_rate[1, ])
  z[1] <- Z[1, ] %*% 1:s
  
  #2期目以降の潜在変数zを逐次的に割り当てる
  for(i in 2:n){
    #多項分布より潜在変数zを生成
    index <- index_s[[z[i-1]]]
    z_rate[i, index] <- Pr[z[i-1], index]*LLi[i, index] / sum(Pr[z[i-1], index]*LLi[i, index])
    Z[i, ] <- extraDistr::rmnom(1, 1, z_rate[i, ])
    z[i] <- Z[i, ] %*% 1:s
  }
  
  ##マルコフ推移行列のパラメータをサンプリング
  index_z <- list()
  for(j in 1:s){
    index_z[[j]] <- which(z[1:(n-1)]==j) + 1  
    index <- index_z[[j]]   #1期前の潜在変数の割当を特定
    par <- colSums(Z[index, index_s[[j]]])   #ベータ分布のパラメータ
    pr <- rbeta(1, par[1]+phi01[1], par[2]+phi01[2])   #ベータ分布よりパラメータをサンプリング
    Pr[j, index_s[[j]]] <- c(pr, 1-pr)
  }
  r <- as.numeric(extraDistr::rdirichlet(1, colSums(Z) + phi02))   #ディクレリ分布より混合率をサンプリング
  
  
  ##観測モデルのパラメータをサンプリング
  mu <- rep(0, n)   #回帰モデルの平均構造の格納用配列
  
  for(j in 1:s){
    #潜在変数の割当を抽出
    index <- which(z==j)
    
    #線形回帰モデルの回帰係数をサンプリング
    x <- X[index, ]
    XXV <- solve(t(x) %*% x + cov0)
    Xy <- t(x) %*% y[index, 1]
    beta_mu <- XXV %*% (Xy + cov0 %*% alpha0)
    oldbeta[, j] <- mvrnorm(1, beta_mu, oldsigma*XXV)   #多変量正規分布より回帰係数をサンプリング
    mu[index] <- x %*% oldbeta[, j]
    
    #ベータ分布のパラメータをサンプリング
    taud <- oldtau[j, ]
    taun <- abs(taud + mvrnorm(1, rep(0, 2), diag(rw[j, ], 2)))
    
    #対数尤度と対数事前分布を計算
    lognew <- sum(dbeta(y[index, 2], taun[1], taun[2], log=TRUE))
    logold <- sum(dbeta(y[index, 2], taud[1], taud[2], log=TRUE))
    logpnew <- lndMvn(taun, lambda0, omega0)
    logpold <- lndMvn(taud, lambda0, omega0)
    
    #MHサンプリング
    alpha <- min(1, exp(lognew + logpnew - logold - logpold))
    if(alpha == "NAN") alpha <- -1
    
    #一様乱数を発生
    u <- runif(1)
    
    #u < alphaなら新しいパラメータを採択
    if(u < alpha){
      oldtau[j, ] <- taun
      
      #そうでないならパラメータを更新しない
    } else {
      oldtau[j, ] <- taud
    }
  }
  
  #線形回帰モデルの共通の分散をサンプリング
  er <- y[, 1] - mu
  s1 <- s0 + t(er) %*% er
  v1 <- v0 + n
  oldsigma <- sqrt(1/(rgamma(1, v1/2, s1/2)))   #逆ガンマ分布からsigma^2をサンプリング
  
  
  ##サンプリング結果の保存と表示
  if(rp%%keep==0){
  
    #サンプリング結果の格納
    mkeep <- rp/keep
    BETA[mkeep, ] <- as.numeric(oldbeta)
    SIGMA[mkeep] <- oldsigma
    LAMBDA[mkeep, ] <- as.numeric(oldtau)
    Zi[mkeep, ] <- z
    PR[, , mkeep] <- Pr
    
    #サンプリング結果の表示
    print(rp)
    print(sum(log(LLi) * Z))
    print(round(cbind(Pr, Pr0), 3))
    print(round(cbind(oldbeta, beta0), 3))
    print(round(c(oldsigma, sigma0), 3))
    print(round(cbind(oldtau, tau0), 3))
  }
}

####サンプリング結果の要約と可視化####
burnin <- 2000/keep
RS <- R/keep

##サンプリング結果の可視化
matplot(BETA[, 1:3], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA[, 4:6], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA[, 7:9], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(BETA[, 10:12], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
plot(1:RS, SIGMA, type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(LAMBDA[, 1:s], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(LAMBDA[, (s+1):ncol(LAMBDA)], type="l", xlab="サンプリング回数", ylab="パラメータ推定値")
matplot(t(PR[1, , ]), type="l", xlab="サンプリング回数", ylab="推移確率のサンプリング結果")
matplot(t(PR[2, , ]), type="l", xlab="サンプリング回数", ylab="推移確率のサンプリング結果")
matplot(t(PR[3, , ]), type="l", xlab="サンプリング回数", ylab="推移確率のサンプリング結果")
matplot(t(PR[4, , ]), type="l", xlab="サンプリング回数", ylab="推移確率のサンプリング結果")

##サンプリング結果の要約
#線形回帰の回帰係数の推定値
round(cbind(matrix(colMeans(BETA[burnin:RS, ]), nrow=k1+1, ncol=s), beta0), 3)
round(matrix(apply(BETA[burnin:RS, ], 2, sd), nrow=k1+1, ncol=s), 3)

#線形回帰の標準偏差の推定値
round(c(mean(SIGMA[burnin:RS]), sigma0), 3)
round(sd(SIGMA[burnin:RS]), 3)

#ベータ分布のパラメータの推定値
round(cbind(matrix(colMeans(LAMBDA[burnin:RS, ]), nrow=s, ncol=2), tau0), 3)
round(matrix(apply(LAMBDA[burnin:RS, ], 2, sd), nrow=s, ncol=2), 3)

#マルコフ推移確率の推定値
round(cbind(apply(PR[, , burnin:RS], c(1, 2), mean), Pr0), 3)
round(apply(PR[, , burnin:RS], c(1, 2), sd), 3)

##予測精度を確認
W <- t(apply(rbind(Zi, matrix(1:s, nrow=s, ncol=n)), 2, table)) - 1
W_rate <- W / rowSums(W)   #潜在変数の割当確率
w <- apply(W_rate, 1, which.max)   #潜在変数の割当
round(cbind(W_rate, w, z0), 3)   #真の潜在変数との比較
sum(w==z0) / n   #予測の正答率


