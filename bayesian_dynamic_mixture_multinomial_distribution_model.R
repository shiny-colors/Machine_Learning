#####ベイジアン潜在マルコフ推移多項分布モデル#####
library(MASS)
library(flexmix)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(2578)

####データの発生####
##データの設定
k <- 10   #セグメント数
hh <- 5000   #ユーザー数
item <- 500   #アイテム数
pt <- rpois(hh, rgamma(hh, 15, 0.7))   #期間数
pt[pt < 5] <- ceiling(runif(sum(pt < 5), 5, 10))
hhpt <- sum(pt)   #総レコード数
s <- rpois(hhpt, rgamma(hhpt, 11.5, 0.8))   #アイテム購買数
s[s < 3] <- ceiling(runif(sum(s < 3), 3, 10))

#IDの設定
u_id <- rep(1:hh, pt)
t_id <- c()
for(i in 1:hh) {t_id <- c(t_id, 1:pt[i])}
ID <- data.frame(no=1:hhpt, u_id, t_id)


##パラメータを設定
#ディレクリ分布のパラメータ
alpha01 <- seq(3.0, 0.2, length=k*5)[((1:(k*5))%%5)==0]
alpha02 <- matrix(0.3, nrow=k, ncol=k)
diag(alpha02) <- 3.5
alpha11 <- rep(0.3, item)

#ディレクリ分布よりパラメータを生成
omegat <- omega <- extraDistr::rdirichlet(1, alpha01)   #ユーザーの1期目のセグメント
gammat <- gamma <- extraDistr::rdirichlet(k, alpha02)   #マルコフ推移行列
thetat <- theta <- extraDistr::rdirichlet(k, alpha11)   #アイテム購買のパラメータ


##ユーザーごとにアイテム購買行列を生成する
Data <- matrix(0, nrow=hhpt, ncol=item)
Z_list <- list()

for(i in 1:hh){
  if(i%%100==0){
    print(i)
  }
  z_vec <- rep(0, s[i])
  
  for(j in 1:pt[i]){
    
    ##期間ごとにセグメントを生成
    index <- which(u_id==i)[j]
    freq <- s[index]
    
    if(j==1){
      z <- rmnom(1, 1, omega)
      z_vec[j] <- as.numeric(z %*% 1:k)
    } else {
      z <- rmnom(1, 1, gamma[z_vec[j-1], ])
      z_vec[j] <- as.numeric(z %*% 1:k)
    }
    
    ##セグメントに基づきアイテム購買行列を生成
    wn <- colSums(rmnom(freq, 1, theta[z_vec[j], ]))
    Data[index, ] <- wn
    Z_list[[index]] <- z_vec[j]
  }
}

#リスト形式を変換
z <- unlist(Z_list)


####マルコフ連鎖モンテカルロ法で潜在マルコフ推移多項分布モデルを推定####
#対数尤度の目標値
LLst <- sum(dmnom(Data, rowSums(Data), colSums(Data)/sum(Data), log=TRUE))

##アルゴリズムの設定
R <- 10000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##事前分布の設定
#ハイパーパラメータの事前分布
alpha01 <- 1 
alpha02 <- 1
beta01 <- 0.5

##パラメータの初期値
theta <- thetat
r0 <- omegat
r1 <- gammat

tf  <- colSums(Data)/sum(Data)*item
theta <- extraDistr::rdirichlet(k, tf)   #アイテム購買確率の初期値
r0 <- rep(1/k, k)
par <- matrix(0.3, nrow=k, ncol=k)
diag(par) <- 2.5
r1 <- extraDistr::rdirichlet(k, par)


##パラメータの格納用配列
THETA <- array(0, dim=c(k, item, R/keep))
R0 <- matrix(0, nrow=R/keep, ncol=k)
R1 <- array(0, dim=c(k, k, R/keep))
SEG <- matrix(0, nrow=hhpt, ncol=k)
storage.mode(SEG) <- "integer"


##MCMC推定用配列
max_time <- max(t_id)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
for(j in 2:max_time){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
}
Data_const <- lfactorial(s) - rowSums(lfactorial(Data))   #多項分布の密度関数の対数尤度の定数
sparse_data <- as(Data, "CsparseMatrix")


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##レコードごとにセグメントをサンプリング
  #セグメントごとの尤度を推定
  theta_log <- log(t(theta))
  LLi0 <- as.matrix(Data_const + sparse_data %*% theta_log)
  LLi_max <- apply(LLi0, 1, max)
  LLi <- exp(LLi0 - LLi_max)
  
  #セグメント割当確率の推定とセグメントの生成
  z_rate <- matrix(0, nrow=hhpt, ncol=k)
  Zi <- matrix(0, nrow=hhpt, ncol=k)
  z_vec <- rep(0, hhpt)
  rf02 <- matrix(0, nrow=k, ncol=k) 
  
  for(j in 1:max_time){
    if(j==1){
      #セグメントの割当確率
      LLs <- matrix(r0, nrow=length(index_t11), ncol=k, byrow=T) * LLi[index_t11, ]   #重み付き尤度
      z_rate[index_t11, ] <- LLs / rowSums(LLs)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi[index_t11, ] <- rmnom(length(index_t11), 1, z_rate[index_t11, ])
      z_vec[index_t11] <- as.numeric(Zi[index_t11, ] %*% 1:k)
      
      #混合率のパラメータを更新
      rf01 <- colSums(Zi[index_t11, ])
      
    } else {
      
      #セグメントの割当確率
      index <- index_t22[[j]]
      z_vec[index_t21[[j]]]
      LLs <- r1[z_vec[index_t21[[j]]], , drop=FALSE] * LLi[index, , drop=FALSE]   #重み付き尤度
      z_rate[index, ] <- LLs / rowSums(LLs)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi[index, ] <- rmnom(length(index), 1, z_rate[index, ])
      z_vec[index] <- as.numeric(Zi[index, ] %*% 1:k)
      
      #混合率のパラメータを更新
      rf02 <- rf02 + t(Zi[index_t21[[j]], , drop=FALSE]) %*% Zi[index, , drop=FALSE]   #マルコフ推移
    }
  }
 
  #ディクレリ分布から混合率をサンプリング
  rf11 <- colSums(Zi[index_t11, ]) + alpha01
  rf12 <- rf02 + alpha01
  r0 <- extraDistr::rdirichlet(1, rf11)
  r1 <- extraDistr::rdirichlet(k, rf12)
  
  #単語分布psiをサンプリング
  wf0 <- matrix(0, nrow=k, ncol=item)
  for(j in 1:k){
    wf0[j, ] <- colSums(sparse_data * Zi[, j])
  }
  wf <- wf0 + alpha01
  theta <- extraDistr::rdirichlet(k, wf)

  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    THETA[, , mkeep] <- theta
    R0[mkeep, ] <- r0
    R1[, , mkeep] <- r1
    
    #トピック割当はバーンイン期間を超えたら格納する
    if(mkeep >= burnin & rp%%keep==0){
      SEG <- SEG + Zi
    }
    
    #サンプリング結果を確認
    if(rp%%disp==0){
      print(rp)
      print(round(cbind(theta[, 1:10], thetat[, 1:10]), 3))
      print(round(cbind(r1, gamma), 3))
    }
  }
}

####サンプリング結果の可視化と要約####
burnin <- 1000/keep   #バーンイン期間
RS <- R/keep

##サンプリング結果の可視化
#アイテム購買確率のサンプリング結果
matplot(t(THETA[, 1, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 50, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 100, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 150, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 200, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 250, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 300, ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(THETA[, 350, ]), type="l", xlab="サンプリング数", ylab="パラメータ")

#マルコフ推移確率のサンプリング結果の可視化
matplot(R0, type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[1, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[2, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[3, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[4, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[5, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[6, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[7, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[8, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[9, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")
matplot(t(R1[10, , ]), type="l", xlab="サンプリング数", ylab="パラメータ")


##サンプリング結果の要約推定量
#セグメント分布の事後推定量
seg_mu <- SEG / rowSums(SEG)
segment <- apply(seg_mu, 1, which.max)
round(cbind(z, seg=segment, seg_mu), 3)   #セグメント割当と真のセグメントの比較

#マルコフ推移確率の事後推定量
round(rbind(colMeans(R0[burnin:RS, ]), omegat), 3)   #1期目の混合率の事後平均
round(cbind(apply(R1[, , burnin:RS], c(1, 2), mean), gammat), 3)   #マルコフ推移確率の事後平均

#アイテム確率の事後推定量
item_mu <- apply(THETA[, , burnin:RS], c(1, 2), mean)   #アイテム購買確率の事後平均
round(cbind(t(item_mu), t(thetat)), 3)

##対数尤度の比較
LLi <- sum(Data_const + rowSums(seg_mu * Data %*% log(t(item_mu))))   #潜在推移多項分布モデルの周辺対数尤度
c(LLi, LLst)



