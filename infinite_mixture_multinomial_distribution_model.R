#####無限混合多項分布モデル#####
library(MASS)
library(Matrix)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(1853)

####データの発生####
#データの設定
N <- 10000   #サンプル数
seg <- 10   #セグメント数
k <- 200   #変数数
w <- rpois(N, rgamma(N, 40, 0.8))   #サンプルあたりの頻度
hist(w, col="grey", breaks=25, main="アイテム購買数の分布", xlab="アイテム購買数")


#セグメントの設定
r <- as.numeric(extraDistr::rdirichlet(1, rep(5.0, seg)))
Z <- rmnom(N, 1, r)
seg_id <- as.numeric(Z %*% 1:seg)


####応答変数の発生####
##セグメントごとにパラメータの設定
alpha <- rep(0.2, k)
theta <- thetat <- extraDistr::rdirichlet(seg, alpha)

##多項分布よりデータを発生
Data <- matrix(0, nrow=N, ncol=k)
for(i in 1:seg){
  index <- which(seg_id==i)
  Data[index, ] <- rmnom(length(index), w[index], theta[i, ])
}
colnames(Data) <- 1:k
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")

####マルコフ連鎖モンテカルロ法で無限次元混合多項分布モデルを推定####
##アルゴリズムの設定
R <- 10000
burnin <- 1000
keep <- 2
disp <- 20
sbeta <- 1.5
iter <- 0

##事前分布の設定
tau <- 1   #ディクレリ分布の事前分布
alpha <- 1   #CRPの事前分布
beta <- 1/k

##初期値の設定
seg0 <- 2   #初期セグメントは2つ
r <- c(0.5, 0.5)   #混合率の初期値
par <- rep(0.1, k)   #CRP用のパラメータ

#初期セグメントを設定
z <- matrix(0, nrow=N, ncol=seg0)
out <- kmeans(Data, seg0)   #kmeans法

#セグメント割当
z0 <- out$cluster
for(i in 1:seg0){z[z0==i, i] <- 1}   

#セグメントごとのパラメータ
oldpar0 <- extraDistr::rdirichlet(seg0, par)
oldpar <- (oldpar0+beta) / matrix(rowSums(oldpar0+beta), nrow=seg0, ncol=k, byrow=T)


#パラメータの格納用配列
max_seg <- 20
Zi <- matrix(0, nrow=N, ncol=max_seg)
THETA <- array(0, dim=c(max_seg, k, R/keep))
storage.mode(Z) <- "integer"

#データの設定
const <- lfactorial(rowSums(Data)) - rowSums(lfactorial(Data))


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  ##多項分布の混合尤度を計算
  #パラメータごとに対数尤度を計算
  LLind0 <- const + Data %*% t(log(oldpar))
  
  #新しい潜在変数の尤度の計算と尤度の結合
  par_mean0 <- as.numeric(extraDistr::rdirichlet(1, Data[sample(1:N, 1), ] + beta))
  par_mean <- (par_mean0+beta) / sum(par_mean0+beta)
  LL_new <- const + Data %*% log(par_mean)   #新しい潜在変数の対数尤度
  
  LLi0 <- cbind(LLind0, LL_new)
  LLi <- exp(LLi0 - rowMaxs(LLi0))   #尤度に変換

  ##CRPの計算
  gamma0 <- cbind(matrix(colSums(z), nrow=N, ncol=ncol(z), byrow=T) - z, alpha)
  gamma1 <- LLi * gamma0/(N-1-alpha)
  
  ##多項分布より潜在変数をサンプリング
  z_rate <- gamma1 / rowSums(gamma1)   #潜在変数zの割当確率
  z <- rmnom(N, 1, z_rate)
  z <- z[, colSums(z) > 0]
  
  ##多項分布のパラメータを更新
  dir_par <- t(t(Data) %*% z) + tau   #ディクレリ分布のパラメータ
  oldpar <- extraDistr::rdirichlet(ncol(z), dir_par)   #多項分布からパラメータを生成

  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[1:nrow(oldpar), , mkeep] <- oldpar
    
    #繰り返し数がバーンインを超えたらパラメータを格納
    if(rp >= burnin){
      Zi[, 1:ncol(z)] <- Zi[, 1:ncol(z)] + z
    }   
    
    if(rp%%disp==0){
      print(rp)
      print(colSums(z))
    }
  }
}

####サンプリング結果の可視化と要約####
#バーンイン期間
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##サンプリング結果をプロット
matplot(t(THETA[1, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(THETA[2, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(THETA[3, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(THETA[4, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(THETA[5, , 1:(R/keep)]), type="l", ylab="パラメータ")
matplot(t(THETA[6, , 1:(R/keep)]), type="l", ylab="パラメータ")

##サンプリング結果の事後平均
mcmc_seg <- sum(colSums(Zi) >= N)   #推定されたセグメント数

#潜在変数zの推定量
round(Z_mu <- (Zi/rowSums(Zi))[, colSums(Zi) > 0], 3)   #潜在変数の割当確率
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #混合率

#多項分布のパラメータの推定量
theta_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  theta_mu[i, ] <- colMeans(t(THETA[i, , burnin1:(R/keep)]))
}
round(cbind(t(theta_mu), t(thetat)), 3) #真のパラメータと比較

