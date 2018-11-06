#####潜在推移多項分布モデル#####
library(MASS)
library(mclust)
library(flexmix)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)
library(knitr)

#set.seed(573298)

####データの発生####
n <- 3000   #サンプル数
seg <- 5   #セグメント数
k <- 25   #観測変数数
period <- 5   #観測期間
N <- n*period   #総サンプル数
freq <- rpois(n*seg, 35)   #観測数


##IDの設定
id <- rep(1:n, rep(period, n))
time <- rep(1:period, n)
ID <- data.frame(no=1:(n*period), id=id, time=time)

##セグメントの設定
#セグメント推移確率行列の設定
P_seg0 <- matrix(0, nrow=seg, ncol=seg)
for(i in 1:seg){
  rand <- runif(seg, 0.1, 3.5)
  P_seg0[i, ] <- rand
}
diag(P_seg0) <- runif(seg, 3.5, 25.0)   #対角行列を置き換え
P_seg <- P_seg0 / matrix(rowSums(P_seg0), nrow=seg, ncol=seg)   #確率に置き換え

#セグメントを発生させる
#期間1のセグメントの設定
seg_m <- matrix(0, nrow=n, ncol=seg)
seg_m[, 1] <- rep(1:seg, rep(n/seg, seg))

#期間1～5まで逐次的にセグメントを発生させる
for(j in 2:seg){
  for(i in 1:n){
    seg_m[i, j] <- t(rmultinom(1, 1, P_seg[seg_m[i, j-1], ])) %*% 1:seg
  }
}
seg_vec <- as.numeric(t(seg_m))   #セグメントをベクトルに変更
table(seg_vec)
r_rate0 <- do.call(rbind, tapply(seg_vec, ID$time, function(x) table(x)/sum(table(x))))   #混合率
r_rate0 <- matrix(as.numeric(r_rate0), nrow=period, ncol=seg)


####セグメント割当に基づき応答変数を発生####
##セグメント別の確率を発生
P <- matrix(0, nrow=seg, ncol=k)
for(i in 1:seg){
  p <- rgamma(k, 0.8, 0.2)
  P[i, ] <- p / sum(p)
}

##多項分布から応答変数を発生
Y <- t(apply(cbind(freq, seg_vec), 1, function(x) rmultinom(1, x[1], P[x[2], ])))

#データの確認と集計
round(colSums(Y)/sum(Y), 3)
by(Y, seg_vec, function(x) round(colMeans(x), 2))
by(Y, seg_vec, function(x) summary(x))
by(Y, seg_vec, function(x) round(colSums(x)/sum(x), 2))


####EMアルゴリズムで潜在推移混合多項分布モデルを推定####
##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta, n, Y, r1, r2, seg){
  LLind <- matrix(0, nrow=nrow(Y), ncol=seg)
  for(i in 1:seg){
    Li <- apply(cbind(n, Y), 1, function(x) dmultinom(x[-1], x[1], theta[i, ]))   #多項分布の尤度
    LLind[, i] <- Li
  }
  
  #尤度が桁落ちしていたら、微小な数を足す
  LLind <- ifelse(matrix(apply(LLind, 1, min), nrow=nrow(Y), ncol=seg)==0, LLind+10^-305, LLind)
  
  LLho <- r1 * LLind   #観測データの尤度
  z <- LLho / matrix(apply(LLho, 1, sum), nrow=nrow(Y), ncol=seg)   #潜在変数zの計算
  LLosum <- sum(log(apply(r1 * LLind, 1, sum)))   #観測データの対数尤度の計算
  rval <- list(LLob=LLosum, z=z, LL=LLind)
  return(rval)
}

##インデックスを作成
index_time <- list()
for(j in 1:seg){
  index_time[[j]] <- subset(1:nrow(ID), ID$time==j)
}


##初期値の設定
#確率の初期値
theta <- matrix(0, nrow=seg, ncol=k)
for(j in 1:seg){
  p0 <- colSums(Y)/sum(Y) + runif(k, 0, 0.75)
  theta[j, ] <- p0 / sum(p0)
}

#混合率の初期値
r <- list()
r[[1]] <- c(0.15, 0.2, 0.2, 0.2, 0.25)
R <- r[[1]]

##対数尤度の初期化
y_seg <- Y[index_time[[j]], ]
freq_seg <- freq[index_time[[j]]]
r1_m <- matrix(r[[1]], nrow=n, ncol=seg, byrow=T)

#潜在変数zと対数尤度の計算
Z <- matrix(0, nrow=N, ncol=seg)
Z_list <- list()
L <- LLobz(theta, freq_seg, y_seg, r1_m, r1_m, seg)
z <- L$z
Z[index_time[[1]], ] <- z
LL <- L$LLob
LLsum <- LL*seg*2

##完全データの対数尤度からパラメータを最大化する
for(j in 1:seg){
  thetaseg <- colSums(matrix(z[, j], nrow=nrow(y_seg), ncol=k) * y_seg) / sum(z[, j] * freq_seg)   #重み付き多項分布の最尤推定
  theta[j, ] <-thetaseg
}

#EMアルゴリズムの更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1 
iter <- 0


####EMアルゴリズムで対数尤度を最大化####
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す

  #Eステップで観測データの対数尤度と潜在変数zの計算
  LL <- c()
  for(j in 1:seg){
    y_seg <- Y[index_time[[j]], ]
    freq_seg <- freq[index_time[[j]]]
    
    #期間が1なら1回前の反復の混合率を期間が2～なら1期前の混合率を事前分布とする
    if(j==1){
      r1 <- r[[j]]
      r1_m <- matrix(r1, nrow=n, ncol=seg, byrow=T)
    } else {
      r1 <- r[[j-1]]
      r1_m <- matrix(r1, nrow=n, ncol=seg, byrow=T)
    }
    
    #観測データの対数尤度と潜在変数zの計算
    L <- LLobz(theta, freq_seg, y_seg, r1_m, r1_m, seg)
    Z[index_time[[j]], ] <- L$z
    Z_list[[j]] <- L$z
    LL <- c(LL, L$LLob)
    
    #混合率を計算
    if(j > 1){
      r[[j]] <- colSums(L$z)/n   #混合率の計算
    }
  }
  
  #対数尤度の和と全体の混合率を計算
  LLsum1 <- sum(LL)   
  r[[1]] <- colSums(Z[index_time[[1]], ])/n   #1期目の混合率
  
  ##Mステップで完全データの対数尤度からパラメータを推定
  for(j in 1:seg){
    thetaseg <- colSums(matrix(Z[, j], nrow=nrow(Y), ncol=k) * Y) / sum(Z[, j] * freq)   #重み付き多項分布の最尤推定
    theta[j, ] <-thetaseg
  }
  
  ##EMアルゴリズムのパラメータを更新
  iter <- iter+1   
  dl <- LLsum-LLsum1
  LLsum <- LLsum1
  print(LLsum1)
}

####推定結果と適合度の確認####
##パラメータの計算
round(r_rate <- do.call(rbind, r), 3)   #混合率の推定値
z_vec　<- apply(Z, 1, which.max)   #セグメント割当


##推定されたパラメータと真のパラメータの比較
round(rbind(theta, P), 3)   #多項分布のパラメータ
round(cbind(z_vec, seg_vec))   #セグメント割当
round(cbind(r_rate, r_rate0), 3)   #混合率


##適合度の計算
round(LLsum1, 3)   #最大化された対数尤度
round(AIC <- -2*LLsum1 + 2*(length(theta) + seg), 3)   #AIC
round(BIC <- -2*LLsum1 + log(N)*(length(theta) + seg), 3)   #BIC
