#####Hidden Marcov Stochastic Block Multinomial model#####
library(MASS)
library(Matrix)
library(bayesm)
library(MCMCpack)
library(gtools)
library(extraDistr)
library(matrixStats)
library(reshape2)
library(qrmtools)
library(slfm)
library(HMM)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####データの発生####
#データの設定
hh <- 3000   #ユーザー数
pt <- rpois(hh, 20)   #観測期間
max_pt <- max(pt)
hhpt <- sum(pt)   #総レコード数
item <- 2000   #アイテム数
seg_u <- 8   #ユーザーのセグメント数
seg_i <- 7   #アイテムのセグメント数

#IDを設定
u_id <- rep(1:hh, pt)
t_id <- c()
for(i in 1:hh){
  t_id <- c(t_id, 1:pt[i])
}

#インデックスを作成
u_index <- list()
for(i in 1:hh){u_index[[i]] <- which(u_id==i)}

##パラメータとセグメントを生成
#パラメータを生成
alpha01 <- rep(10, seg_u)
alpha02 <- matrix(2.0, nrow=seg_u, ncol=seg_u)
diag(alpha02) <- runif(seg_u, 10, 16)
gamma01 <- as.numeric(extraDistr::rdirichlet(1, alpha01))
gamma02 <- extraDistr::rdirichlet(seg_u, alpha02)

##ユーザーセグメントを生成
z01 <- matrix(0, nrow=hhpt, ncol=seg_u)
z01_vec <- rep(0, hhpt)
for(i in 1:hh){
  for(j in 1:pt[i]){
    if(j==1){
      #1期目のセグメントを生成
      z01[u_index[[i]][j], ] <- rmnom(1, 1, gamma01)
      z01_vec[u_index[[i]][j]] <- as.numeric(z01[u_index[[i]][j], ] %*% 1:seg_u)
      
    } else {
      
      #2期目以降のセグメントを生成
      z01[u_index[[i]][j], ] <- rmnom(1, 1, gamma02[z01_vec[u_index[[i]][j-1]], ])
      z01_vec[u_index[[i]][j]] <- as.numeric(z01[u_index[[i]][j], ] %*% 1:seg_u)
    }
  }
}
z1_cnt <- colSums(z01)

##アイテムセグメントを生成
alpha03 <- rep(20, seg_i)
omega01 <- as.numeric(extraDistr::rdirichlet(1, alpha03))
z02 <- rmnom(item, 1, omega01)
z02_vec <- as.numeric(z02 %*% 1:seg_i)
z2_cnt <- colSums(z02)


##観測モデルのパラメータを生成
#ユーザーセグメント×アイテムセグメントのベータ事前分布のパラメータを発生
hist(rbeta(10000, 0.4, 4.5), col="grey", breaks=25, main="ベータ分布からの乱数", xlab="パラメータ")
theta0 <- matrix(rbeta(seg_u*seg_i, 0.4, 4.5), nrow=seg_u, ncol=seg_i)
round(theta0, 3)

##ベルヌーイ分布から観測行列を発生させる
#ベルヌーイ分布から共起行列を生成
Data <- matrix(0, nrow=hhpt, ncol=item)

for(i in 1:seg_u){
  print(i)
  for(j in 1:seg_i){
    n <- z1_cnt[i] * z2_cnt[j]
    z1_cnt[i]
    Data[z01_vec==i, z02_vec==j] <- matrix(rbinom(n, 1, theta0[i, j]), nrow=z1_cnt[i], ncol=z2_cnt[j])
  }
}

Data_T <- t(Data)
storage.mode(Data) <- "integer"
storage.mode(Data_T) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")
sparse_data_T <- as(Data_T, "CsparseMatrix")  
gc(); gc()

#多項分布のパラメータを設定
thetat <- matrix(0, nrow=seg_u, ncol=seg_i)
for(i in 1:seg_u){
  n <- sum(sparse_data[z01_vec==i, ])
  for(j in 1:seg_i){
    thetat[i, j] <- sum(sparse_data[z01_vec==i, z02_vec==j]) / n
  }
}

phit <- matrix(0, nrow=seg_i, nco=seg_u)
for(i in 1:seg_i){
  n <- sum(sparse_data[, z02_vec==i])
  for(j in 1:seg_u){
    phit[i, j] <- sum(sparse_data[z01_vec==j, z02_vec==i]) / n
  }
}


####マルコフ連鎖モンテカルロ法でHM Sparse Stochastic Block modelを推定####
##アルゴリズムの設定
R <- 5000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##事前分布の設定
tau1 <- rep(1, seg_u)
tau2 <- matrix(1, nrow=seg_u, ncol=seg_u)
diag(tau2) <- 5
tau3 <- rep(1, seg_i)
alpha1 <- matrix(1, nrow=seg_u, ncol=seg_i)   #ユーザーセグメントのディクレリ事前分布
alpha2 <- matrix(1, nrow=seg_i, ncol=seg_u)   #アイテムセグメントのディクレリ事前分布


##初期値の設定
#混合率の初期値
alpha01 <- rep(20, seg_u)
alpha02 <- matrix(4, nrow=seg_u, ncol=seg_u)
alpha03 <- rep(20, seg_i)
diag(alpha02) <- 20
gamma1 <- extraDistr::rdirichlet(1, alpha01)
gamma2 <- extraDistr::rdirichlet(seg_u, alpha02)
omega1 <- extraDistr::rdirichlet(1, alpha03)

#ブロックごとのパラメータの初期値
index_u <- floor(seq(1, hhpt, length=seg_u+1))
index_i <- floor(seq(1, item, length=seg_i+1))
sortlist1 <- order(rowSums(sparse_data))
sortlist2 <- order(colSums(sparse_data), decreasing=TRUE)

#クラス割当の初期値を設定
z1 <- rep(0, hhpt)
z2 <- rep(0, item)

for(i in 1:(length(index_u)-1)){
  #ユーザーのクラス割当のインデックスを設定
  index1 <- sortlist1[index_u[i]:index_u[i+1]]
  z1[index1] <- i
  
  for(j in 1:(length(index_i)-1)){
    #アイテムのクラス割当のインデックスを設定
    index2 <- sortlist2[index_i[j]:index_i[j+1]]
    z2[index2] <- j
  }
}

#セグメントインデックスを作成
index1 <- list()
index2 <- list()
user_vec <- matrix(0, nrow=hhpt, ncol=seg_u)
item_vec <- matrix(0, nrow=item, ncol=seg_i)
n1 <- c()
n2 <- c()

for(i in 1:seg_u){
  index1[[i]] <- which(z1==i)
  user_vec[index1[[i]], i] <- 1
  n1 <- c(n1, sum(sparse_data[index1[[i]], ]))
}
for(j in 1:seg_i){
  index2[[j]] <- which(z2==j)
  item_vec[index2[[j]], j] <- 1
  n2 <- c(n2, sum(sparse_data[, index2[[j]]]))
}

#パラメータの初期値を設定
#アイテムセグメントの初期パラメータ
oldtheta <- matrix(0, nrow=seg_u, ncol=seg_i) 
for(i in 1:seg_u){
  for(j in 1:(seg_i-1)){
    freq <- sum(sparse_data[index1[[i]], index2[[j]]])
    oldtheta[i, j] <- freq / n1[i]
  }
}
oldtheta[, seg_i] <- 1 - rowSums(oldtheta)
oldtheta <- (oldtheta + 0.00001) / rowSums(oldtheta + 0.00001)

#アイテムセグメントの初期パラメータ
oldphi <- matrix(0, nrow=seg_i, ncol=seg_u)
for(i in 1:seg_i){
  for(j in 1:(seg_u-1)){
    freq <- sum(sparse_data[index1[[j]], index2[[i]]])
    oldphi[i, j] <- freq / n2[i]
  }
}
oldphi[, seg_u] <- 1 - rowSums(oldphi)
oldphi <- (oldphi + 0.00001) / rowSums(oldphi + 0.00001)

#ユーザー、アイテムごとの購買数
n_user <- rowSums(sparse_data)
n_item <- colSums(sparse_data)

##パラメータの格納用配列
THETA <- array(0, dim=c(seg_u, seg_i, R/keep))
PHI <- array(0, dim=c(seg_i, seg_u, R/keep))
GAMMA1 <- matrix(0, nrow=R/keep, ncol=seg_u)
GAMMA2 <- array(0, dim=c(seg_u, seg_u, R/keep))
OMEGA1 <- matrix(0, nrow=R/keep, ncol=seg_i)
SEG1 <- matrix(0, hhpt, ncol=seg_u)
SEG2 <- matrix(0, item, ncol=seg_i)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()


##インデックスを作成
max_pt <- max(pt)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
n_time <- rep(0, max_pt)
n_time[1] <- length(index_t11)
for(j in 2:max_pt){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
  n_time[j] <- length(index_t22[[j]])
}


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){

  ##ユーザーごとのセグメント割当を生成
  #ユーザーのセグメント割当ごとの尤度を推定
  y1 <- as.matrix(sparse_data %*% item_vec)   #アイテムセグメントごとの購入頻度
  LLi0 <- y1 %*% t(log(oldtheta))   #セグメントごとの多項分布の対数尤度
  LLi_max <- rowMaxs(LLi0)
  LLi1 <- exp(LLi0 - LLi_max)   #尤度に変換
  
  #セグメント割当確率の推定とセグメントの生成
  z_rate1 <- matrix(0, nrow=hhpt, ncol=seg_u)
  Zi1 <- matrix(0, nrow=hhpt, ncol=seg_u)
  z1_vec <- rep(0, hhpt)
  rf02 <- matrix(0, nrow=seg_u, ncol=seg_u)
  
  for(j in 1:max_pt){
    if(j==1){
      #セグメント割当確率を推定
      n <- n_time[j]
      LLs <- matrix(gamma1, nrow=n, ncol=seg_u, byrow=T) * LLi1[index_t11, ]   #重み付き尤度
      matrix(gamma1, nrow=n, ncol=seg_u, byrow=T)
      z_rate1[index_t11, ] <- LLs / rowSums(LLs)   #割当確率
      
      #多項分布からセグメントを生成
      Zi1[index_t11, ] <- rmnom(n, 1, z_rate1[index_t11, ])
      z1_vec[index_t11] <- as.numeric(Zi1[index_t11, ] %*% 1:seg_u)
      
      #混合率のパラメータを更新
      rf01 <- colSums(Zi1[index_t11, ])
      
    } else {
      
      #セグメントの割当確率
      index <- index_t22[[j]]
      n <- n_time[j]
      LLs <- gamma2[z1_vec[index_t21[[j]]], , drop=FALSE] * LLi1[index, , drop=FALSE]   #重み付き尤度
      z_rate1[index, ] <- LLs / rowSums(LLs)   #割当確率
      
      #多項分布よりセグメントを生成
      Zi1[index, ] <- rmnom(n, 1, z_rate1[index, ])
      z1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:seg_u)
      
      #混合率のパラメータを更新
      rf02 <- rf02 + t(Zi1[index_t21[[j]], , drop=FALSE]) %*% Zi1[index, , drop=FALSE]   #マルコフ推移
    }
  }
  user_vec <- Zi1
  
  #マルコフ推移行列のパラメータを更新
  rf1 <- rf01 + tau1
  rf2 <- rf02 + tau2
  gamma1 <- as.numeric(extraDistr::rdirichlet(1, rf1))
  gamma2 <- extraDistr::rdirichlet(seg_u, rf2)
  
  
  ##アイテムごとにセグメント割当を生成
  #セグメントごとに多項分布の尤度を計算
  y2 <- sparse_data_T %*% user_vec   #ユーザセグメントごとの購買頻度
  LLi2 <- as.matrix(y2 %*% log(t(oldphi)))
  
  #logsumexpの尤度を計算
  LLi_max <- rowMaxs(LLi2)
  r_matrix <- matrix(omega1, nrow=item, ncol=seg_i, byrow=T)
  
  #セグメント割当確率の推定とセグメントの生成
  expl <- r_matrix * exp(LLi2 - LLi_max)
  z2_rate <- expl / rowSums(expl)   #セグメント割当確率
  item_vec <- Zi2 <- rmnom(item, 1, z2_rate)   #多項分布からセグメント割当を生成
  z2_vec <- as.numeric(Zi2 %*% 1:seg_i)

  #混合率の更新
  rf3 <- colSums(Zi2) + tau3
  omega1 <- as.numeric(extraDistr::rdirichlet(1, rf3))
  
  
  ##セグメントインデックスを作成
  index1 <- list()
  index2 <- list()
  for(i in 1:seg_u){index1[[i]] <- which(z1_vec==i)}
  for(j in 1:seg_i){index2[[j]] <- which(z2_vec==j)}

  
  ##ユーザーおよびアイテムのパラメータをサンプリング
  #ディクレリ分布からユーザーセグメントのパラメータをサンプリング
  freq_user <- matrix(0, nrow=seg_u, ncol=seg_i)
  for(i in 1:seg_u){
    x <- sparse_data[index1[[i]], , drop=FALSE]
    for(j in 1:seg_i){
      freq_user[i, j] <- sum(x[, index2[[j]]])
    }
  }
  freq_item <- t(freq_user)   #アイテムのパラメータはユーザーパラメータを反転させるだけ

  
  #ディクレリ分布からアイテムセグメントのパラメータをサンプリング
  oldtheta <- extraDistr::rdirichlet(seg_u, freq_user + alpha1)
  oldphi <- extraDistr::rdirichlet(seg_i, freq_item + alpha2)
  
  ##サンプリング結果の格納と表示
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta 
    PHI[, , mkeep] <- oldphi
    GAMMA1[mkeep, ] <- gamma1
    GAMMA2[, , mkeep] <- gamma2
    OMEGA1[mkeep, ] <- omega1
    
    if(rp >= 500){
    SEG1 <- SEG1 + Zi1 
    SEG2 <- SEG2 + Zi2 
    }
    
    if(rp%%disp==0){
      print(rp)
      print(round(rbind(gamma1, gamma01), 3))
      print(round(cbind(gamma2, gamma02), 3))
      print(round(cbind(oldtheta, thetat), 3))
      print(round(cbind(oldphi, phit), 3))
    }
  }
}

