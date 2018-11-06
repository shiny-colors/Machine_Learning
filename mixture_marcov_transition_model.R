#####有限混合マルコフ推移モデル####
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
hh <- 5000   #ユーザー数
S <- 8   #ページ数
k <- 5   #混合数
seg <- as.numeric(rmnom(hh, 1, rep(1, k)) %*% 1:k)

##パラメータの設定
#マルコフ推移行列の設定
Pr <- array(0, dim=c(S-1, S, k))
for(i in 1:k){
  for(j in 1:(S-1)){
    if(j==1){
      Pr[j, -c(j, S), i] <- extraDistr::rdirichlet(1, rep(1, S-2))   
    } else {
      Pr[j, -j, i] <- extraDistr::rdirichlet(1, rep(1, S-1))
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
r <- rep(0.2, k)

####EMアルゴリズムで有限混合マルコフ推移モデルを推定####
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

##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta, r, Data, id, index, hh, k){

  #潜在変数ごとの尤度を計算
  LLind <- matrix(0, nrow=hh, ncol=k)
  for(j in 1:k){
    Li <- rowProds(theta[index, , j] ^ Data)
    LLind[, j] <- tapply(Li, id, prod)
  }
  
  #潜在変数の割当確率を計算
  LLho <- matrix(r, nrow=hh, ncol=k, byrow=T) * LLind   #観測データの尤度
  z <- LLho / matrix(rowSums(LLho), nrow=hh, ncol=k)   #潜在変数zの割当確率
  LL <- sum(log(rowSums(LLho)))   #観測データの対数尤度
  rval <- list(LLob=LL, z=z, LL=LLind)
  return(rval)
}

##初期値の設定
#パラメータの初期値
theta <- array(0, dim=c(S-1, S, k))
for(i in 1:k){
  for(j in 1:(S-1)){
    if(j==1){
      theta[j, -c(j, S), i] <- extraDistr::rdirichlet(1, rep(5, S-2))   
    } else {
      theta[j, -j, i] <- extraDistr::rdirichlet(1, rep(5, S-1))
    }
  }
}

#混合率の初期値
r <- rep(1/k, k)

#対数尤度の初期化
L <- LLobz(theta, r, Data, id, index_trans, hh, k)
LL1 <- L$LLob
z <- L$z

#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1  

##EMアルゴリズムで有限混合マルコフ推移モデルのパラメータを更新
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #潜在変数zの出力
  z <- L$z   
  
  #Mステップの計算と最適化
  #thetaの推定
  theta <- array(0, dim=c(S-1, S, k))
  
  for(j in 1:k){
    #完全データの対数尤度からthetaの推定量を計算
    wt_data <- matrix(z[id, j], nrow=nrow(Data), ncol=S) * Data   #重み付きデータを作成
    theta0 <- as.matrix(data.frame(id=index_trans, data=wt_data) %>%
                          dplyr::group_by(id) %>%
                          dplyr::summarise_all(funs(sum)))[, 2:(S+1)]
    theta[, , j] <- theta0 / matrix(rowSums(theta0), nrow=S-1, ncol=S)
  }
  
  #混合率を推定
  r <- colSums(z[id, ]) / nrow(Data)
  
  #観測データの対数尤度を計算(Eステップ)
  L <- LLobz(theta, r, Data, id, index_trans, hh, k)
  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####推定結果の確認####
#パラメータの推定値
round(theta, 3)   #マルコフ推移確率の推定値
round(Pr, 3)   #マルコフ推移確率の真値
round(rbind(r1=r, r0=table(seg)/sum(table(seg))), 3)   #混合率

#潜在変数の割当
round(z, 3)   #潜在変数zの割当確率
cbind(seg1=apply(z, 1, which.max), seg0=seg)   #推定された潜在変数と真の潜在変数

#適合度
LL   #最大化された対数尤度
-2*(LL) + 2*(sum(theta[, , 1] > 0)*k) #AIC





