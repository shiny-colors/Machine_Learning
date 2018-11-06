#####ベイジアン無限潜在推移二項分布モデル#####
library(MASS)
library(lda)
library(RMeCab)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(45327)

####データの発生####
##データの設定
hh <- 5000   #ユーザー数
pt <- 7   #観測期間数数
hhpt <- hh*pt   #総サンプル数
seg <- 7   #セグメント数
k <- 50   #観測変数数

##IDの設定
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id, time)

##セグメントの設定
#セグメント推移確率行列の設定
P_seg0 <- matrix(0, nrow=seg, ncol=seg)
for(i in 1:seg){
  rand <- runif(seg, 0.1, 4)
  P_seg0[i, ] <- rand
}
diag(P_seg0) <- runif(seg, 5, 25)   #対角行列を置き換え
P_seg <- P_seg0 / matrix(rowSums(P_seg0), nrow=seg, ncol=seg)   #確率に置き換え

#セグメントを発生させる
#期間1のセグメントを設定
seg_m <- matrix(0, nrow=hh, ncol=pt)
seg_m[, 1] <- rep(1:seg, rep(ceiling(hh/seg), seg))[1:hh]

#期間1~7まで逐次的にセグメントを発生させる
for(j in 2:pt){
  for(i in 1:hh){
    seg_m[i, j] <- extraDistr::rmnom(1, 1, P_seg[seg_m[i, j-1], ]) %*% 1:seg
  }
}
seg_vec <- as.numeric(t(seg_m))   #セグメントをベクトルに変更
r_rate0 <- do.call(rbind, tapply(seg_vec, ID$time, function(x) table(x)/sum(table(x))))   #混合率の真値
r_rate0 <- matrix(as.numeric(r_rate0), nrow=pt, ncol=seg)

####セグメント割当に基づき応答変数を発生####
##セグメント別のパラメータを設定
Pr0 <- matrix(0, nrow=seg, ncol=k)
for(i in 1:seg){
  for(j in 1:k){
    tau1 <- runif(1, 0.5, 4.5)
    tau2 <- runif(1, 2.0, 9.5)
    Pr0[i, j] <- rbeta(1, tau1, tau2)
  }
}


##二項分布から応答変数を発生
Data <- matrix(0, nrow=hhpt, ncol=k)
for(i in 1:seg){
  Data[seg_vec==i, ] <- matrix(rbinom(sum(seg_vec==i)*k, 1, Pr0[i, ]), nrow=sum(seg_vec==i), ncol=k, byrow=T)
}

#データの確認と集計
colSums(Data); colMeans(Data)
by(Data, seg_vec, function(x) round(colMeans(x), 3))
by(Data, seg_vec, function(x) round(colSums(x), 3))
by(Data, seg_vec, function(x) summary(x))


####マルコフ連鎖モンテカルロ法で無限潜在推移確率モデルを推定####
##アルゴリズムの設定
R <- 20000
keep <- 4 
sbeta <- 1.5
iter <- 0

##事前分布の設定
tau <- c(1, 1)   #ベータ分布の事前分布
alpha <- 1   #CRPの事前分布

##初期値の設定
seg0 <- 2   #初期セグメントは2つ
r <- c(0.5, 0.5)   #混合率の初期値
par_mean <- colMeans(Data)   #CRP用のパラメータ

#初期セグメントの設定
z <- matrix(0, nrow=hhpt, ncol=seg0)
z0 <- rbinom(hhpt, 1, 0.5) + 1
for(i in 1:seg0){z[z0==i, i] <- 1}
z_vec <- z %*% 1:seg0

#セグメントごとのパラメータ
oldtheta <- matrix(rbeta(k*seg0, 2.0, 2.0), nrow=seg0, ncol=k)

##パラメータの格納用
max_seg <- 20
Z <- matrix(0, nrow=hhpt, ncol=max_seg)
P <- array(0, dim=c(max_seg, k, R/keep))
storage.mode(Z) <- "integer"


####MCMCでパラメータをサンプリング####
for(rp in 1:R){
  
  LLind0 <- matrix(0, nrow=hhpt, ncol=nrow(oldtheta))
  for(j in 1:ncol(LLind0)){
    Li <- Data %*% log(oldtheta[j, ]) + (1-Data) %*% log(1-oldtheta[j, ])
    LLind0[, j] <- Li
  }
  
  #新しい潜在変数の尤度の計算と尤度の結合
  par <- colMeans(Data)
  LL_new <- Data %*% log(par) + (1-Data) %*% log(1-par)
  LLi0 <- cbind(LLind0, LL_new)
  LLi <- exp(LLi0 - max(LLi0))   #尤度に変換
  
  ##潜在変数の割当確率の事前分布とCRPの計算
  r0 <- list()
  r <- matrix(0, nrow=hhpt, ncol=ncol(z)+1)
  
  #1期では全体の混合率を2期以降は前期の混合率を採用
  for(j in 1:pt){
    if(j==1){
      index <- which(ID$time==j)
      r0[[j]] <- cbind(matrix(colSums(z), nrow=length(index), ncol=ncol(z), byrow=T) - z[index, ], alpha)
      r[index, ] <- r0[[j]] / (hhpt-1-alpha/pt)
    }
    if(j!=1){
      index_now <- which(ID$time==j)
      index_obs <- which(ID$time==j-1)
      r0[[j]] <- cbind(matrix(colSums(z[index_obs, ]), nrow=length(index_obs), ncol=ncol(z), byrow=T) - z[index_obs, ], alpha/pt)
      r[index_now, ] <- r0[[j]] / (length(index_obs)-1-alpha/pt)
    }
  }
  
  ##潜在変数の割当確率の計算と潜在変数のサンプリング
  gamma <- LLi * r
  z_rate <- gamma / rowSums(gamma)   #潜在変数zの割当確率
  z <- rmnom(hhpt, 1, z_rate)   #多項分布から潜在変数zをサンプリング
  z[is.nan(z)] <- 0
  z <- z[, colSums(z) > 2]
  
  ##二項分布のパラメータを更新
  oldtheta <- matrix(0, nrow=ncol(z), ncol=k) 
  for(j in 1:ncol(z)){
    
    #パラメータ更新に必要なデータを抽出
    index <- which(z[, j]==1)
    x <- Data[index, ]   #割当セグメントに該当するデータを抽出
    
    if(length(x) > k){
      y <- colSums(x)
      n <- nrow(x)
    } else {
      y <- x
      n <- 1
    }
    
    #ベータ分布からパラメータを発生
    phi1 <- tau[1] + y
    phi2 <- tau[2] + n - y
    oldtheta[j, ] <- rbeta(k, phi1, phi2)
  }
  
  ##パラメータの格納とサンプリング結果の表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    if(rp >= R/4){Z[, 1:ncol(z)] <- Z[, 1:ncol(z)] + z}   #繰り返し数が最大反復数の1/4を超えたらパラメータを格納
    P[1:nrow(oldtheta), , mkeep] <- oldtheta
    
    print(rp)
    print(round(colMeans(z), 3))
    print(round(rbind(oldtheta, Pr0)[, 1:15], 3))
  }
}


####サンプリング結果の要約と可視化####
burnin <- 2500/keep
RS <- R/keep

##サンプリング結果の可視化
matplot(t(P[1:seg, 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(P[1:seg, 2, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(P[1:seg, 3, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(P[1:seg, 4, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(P[1:seg, 5, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")

##サンプリング結果の要約
round(rbind(apply(P[1:seg, , burnin:RS], c(1, 2), mean), Pr0), 3)   #事後平均
round(apply(P[1:seg, , burnin:RS], c(1, 2), sd), 3)   #事後標準偏差
round(cbind(Z[, colSums(Z) > 0] / rowSums(Z), seg=seg_vec), 3)   #潜在変数の割当の事後分布


