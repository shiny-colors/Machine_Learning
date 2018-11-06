#####ベイジアン有限混合マルコフ推移モデル#####
library(MASS)
library(Matrix)
library(flexmix)
library(mclust)
library(matrixStats)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(90345)

####データの発生####
hh <- 10000   #ユーザー数
S <- 20   #ページ数
k <- 10   #混合数
seg <- as.numeric(rmnom(hh, 1, rep(1, k)) %*% 1:k)
rt <- as.numeric(table(seg)/hh)

##パラメータの設定
#マルコフ推移行列の設定
Pr <- array(0, dim=c(S-1, S, k))
for(i in 1:k){
  for(j in 1:(S-1)){
    if(j==1){
      Pr[j, -S, i] <- extraDistr::rdirichlet(1, rep(1, S-1))   
    } else {
      Pr[j, , i] <- extraDistr::rdirichlet(1, rep(1, S))
    }
  }
}

##ユーザーごとにコンバージョンするまでデータを逐次発生
Data_list <- list()
id_list <- list()

for(i in 1:hh){
  data <- matrix(0, nrow=1000, ncol=S)
  data[1, ] <- rmnom(1, 1, Pr[1, , seg[i]])   #1アクセス目のログを生成
  
  for(j in 2:1000){
    data[j, ] <- rmnom(1, 1, Pr[which(data[j-1, ]==1), , seg[i]])   #2アクセス以降のログを生成
    if(data[j, S]==1) break   #コンバージョンしていたらbreak
  }
  Data_list[[i]] <- data[rowSums(data) > 0, ]
  id_list[[i]] <- rep(i, sum(rowSums(data) > 0))
}
Data <- do.call(rbind, Data_list)
id <- unlist(id_list)
as.matrix(data.frame(id, Data) %>%
            dplyr::group_by(id) %>%
            dplyr::summarize_all(funs(sum)))
r <- rep(1/k, k)


####マルコフ連鎖モンテカルロ法で有限混合マルコフ推移モデルを推定####
##推移ベクトルを作成
index_list <- list()
for(i in 1:hh){
  data <- Data[id==i, ]
  index <- rep(0, nrow(data))
  index[1] <- 1
  index[-1] <- data[1:(nrow(data)-1), ] %*% 1:S
  index_list[[i]] <- index
}
index_trans <- unlist(index_list)

##潜在変数zを計算するための関数
LLobz <- function(theta, r, Data, id, index, hh, k){
  index <- index_trans
  #潜在変数ごとの尤度を計算
  LLind0 <- matrix(0, nrow=hh, ncol=k)
  log_theta <- log(theta)   #パラメータを対数変換
  for(j in 1:k){
    Li <- rowSums(Data * log_theta[index, , j])
    LLind0[, j] <- tapply(Li, id, sum)
  }
  LLind <- exp(LLind0 - apply(LLind0, 1, max))
  
  #潜在変数の割当確率を計算
  LLho <- matrix(r, nrow=hh, ncol=k, byrow=T) * LLind   #観測データの尤度
  z <- LLho / matrix(rowSums(LLho), nrow=hh, ncol=k)   #潜在変数zの割当確率
  rval <- list(z=z)
  return(rval)
}

##アルゴリズムの設定
R <- 2000
keep <- 2  
iter <- 0
burnin <- 200/keep0
disp <- 10

##事前分布の設定
alpha <- 1
beta <- 1

##初期値の設定
#パラメータの初期値
theta <- array(0, dim=c(S-1, S, k))
for(i in 1:k){
  for(j in 1:(S-1)){
    if(j==1){
      theta[j, -S, i] <- extraDistr::rdirichlet(1, rep(5, S-1))   
    } else {
      theta[j, , i] <- extraDistr::rdirichlet(1, rep(5, S))
    }
  }
}
theta <- theta + 0.0001

#混合率の初期値
r <- rep(1/k, k)

##サンプリング結果の保存用配列
THETA <- array(0, dim=c(S-1, S, k, R/keep))
RATE <- matrix(0, nrow=R/keep, ncol=k)
SEG <- matrix(0, nrow=hh, ncol=k)
storage.mode(Z) <- "integer"


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){

  ##ユーザーごとにセグメントをサンプリング
  z <- LLobz(theta, r, Data, id, index_trans, hh, k)$z   #セグメント割当確率
  Zi <- rmnom(hh, 1, z)   #多項分布からセグメントを生成

  ##パラメータをサンプリング
  #混合率をサンプリング
  wsum <- colSums(Zi) + alpha
  r <- as.numeric(extraDistr::rdirichlet(1, wsum))

  #マルコフ推移行列をサンプリング
  theta <- array(0, dim=c(S-1, S, k))
  for(j in 1:k){
    index <- which(Zi[id, j]==1)
    data <- Data[index, ]
    theta0 <- as.matrix(data.frame(id=index_trans[index], data=data) %>%
                          dplyr::group_by(id) %>%
                          dplyr::summarise_all(funs(sum)))[, 2:(S+1)]
    dsum <- theta0 + beta
    theta[, , j] <- extraDistr::rdirichlet(S-1, dsum)
  }
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , , mkeep] <- theta
    RATE[mkeep, ] <- r

    #トピック割当はバーンイン期間を超えたら格納する
    if(mkeep >= burnin & rp%%keep==0){
      SEG <- SEG + Zi
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(round(rbind(r, rt), 3))
    }
  }
}

####推定結果の確認####
burnin <- 200/keep0
RS <- R/keep

##サンプリング結果の可視化
matplot(t(THETA[1, , 1, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[2, , 2, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[3, , 3, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[4, , 4, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(t(THETA[5, , 5, ]), type="l", xlab="サンプリング回数", ylab="パラメータ")
matplot(RATE, type="l", xlab="サンプリング回数", ylab="パラメータ")

#パラメータの推定値
theta_mu <- array(0, dim=c(S-1, S, k))
for(j in 1:k){
  theta_mu[, , j] <- apply(THETA[, , j, burnin:RS], c(1, 2), mean)
}

round(theta_mu, 3)   #マルコフ推移確率の推定値
round(Pr, 3)   #マルコフ推移確率の真値
round(rbind(r1=colMeans(RATE[burnin:RS, ]), rt), 3)   #混合率


#潜在変数の割当
round(Z <- SEG / rowSums(SEG))   #潜在変数zの割当確率
round(cbind(seg0=seg, seg1=apply(Z, 1, which.max), Z), 3)   #推定された潜在変数と真の潜在変数

##適合度
#ユニグラムの対数尤度
LLst1 <- sum(Data %*% log(colMeans(Data)))

#マルコフモデルの対数尤度
par <- matrix(0, nrow=S-1, ncol=S)
for(j in 1:max(index_trans)){
  index <- which(index_trans==j)
  par[j, ] <- colMeans(Data[index, ])
}
log_par <- log(par)
LLi <- Data * log_par[index_trans, ]
LLi[is.nan(LLi)] <- 0
LLst2 <- sum(LLi)

#混合マルコフモデルの対数尤度
LLi <- matrix(0, nrow=nrow(Data), ncol=k)
for(j in 1:k){
  log_theta <- log(theta_mu[, , j])
  LLi[, j] <- rowSums(Z[id, j] * Data * log_theta[index_trans, ])
}
LL <- sum(LLi)

#適合度の比較
round(LLc <- c(LLst1, LLst2, LL), 3)   #対数尤度
round(exp(-LLc / nrow(Data)), 3)   #Perplexity
