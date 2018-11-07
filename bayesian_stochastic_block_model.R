#####ベイジアン確率的ブロックモデル#####
library(MASS)
library(bayesm)
library(MCMCpack)
library(gtools)
library(extraDistr)
library(reshape2)
library(matrixStats)
library(qrmtools)
library(slfm)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####データの発生####
#データの設定
N <- 5000   #ユーザー数
K <- 2000   #アイテム数
seg_u <- 8   #ユーザーのセグメント数
seg_i <- 6   #アイテムのセグメント数

##パラメータとセグメントを発生させる
#ユーザーセグメントを発生
alpha01 <- rep(100, seg_u)
pi01 <- extraDistr::rdirichlet(1, alpha01)
z01 <- t(rmultinom(N, 1, pi01))
z01_vec <- as.numeric(z01 %*% 1:seg_u)
z1_cnt <- as.numeric(table(z01_vec))
mix1 <- colMeans(z01)   #混合率の真値


#アイテムセグメントの発生
alpha02 <- rep(100, seg_i)
pi02 <- extraDistr::rdirichlet(1, alpha02)
z02 <- t(rmultinom(K, 1, pi02))
z02_vec <- as.numeric(z02 %*% 1:seg_i)
z2_cnt <- as.numeric(table(z02_vec))
mix2 <- colMeans(z02)   #混合率の真値


#観測変数のパラメータの設定
#ユーザーセグメント×アイテムセグメントのベータ事前分布のパラメータを発生
hist(rbeta(10000, 1.0, 6.5), col="grey", breaks=25, main="ベータ分布からの乱数", xlab="パラメータ")
theta0 <- matrix(rbeta(seg_u*seg_i, 1.0, 6.5), nrow=seg_u, ncol=seg_i)
round(theta0, 3)


##ベルヌーイ分布から観測行列を発生させる
#ベルヌーイ分布から共起行列を生成
Data <- matrix(0, nrow=N, ncol=K)

for(i in 1:seg_u){
  print(i)
  for(j in 1:seg_i){
    n <- z1_cnt[i] * z2_cnt[j]
    Data[z01_vec==i, z02_vec==j] <- matrix(rbinom(n, 1, theta0[i, j]), nrow=z1_cnt[i], ncol=z2_cnt[j])
  }
}

Data_T <- t(Data)
storage.mode(Data) <- "integer"
storage.mode(Data_T) <- "integer"
sparse_data0 <- as(1-Data, "CsparseMatrix")   
sparse_data1 <- as(Data, "CsparseMatrix")  
sparse_data_T0 <- as(1-Data_T, "CsparseMatrix")   
sparse_data_T1 <- as(Data_T, "CsparseMatrix")   
gc(); gc()


####マルコフ連鎖モンテカルロ法で確率的ブロックモデルを推定####
##アルゴリズムの設定
R <- 5000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##事前分布の設定
alpha1 <- alpha2 <- 1
tau <- c(1, 1)   #ベータ分布の事前分布
alpha1 <- rep(1, seg_u)   #ユーザーセグメントのディクレリ事前分布
alpha2 <- rep(1, seg_i)   #アイテムセグメントのディクレリ事前分布

##初期値の設定
#混合率の初期値
r1 <- rep(1/seg_u, seg_u)
r2 <- rep(1/seg_i, seg_i)

#ブロックごとのパラメータの初期値
index_u <- floor(seq(1, N, length=seg_u+1))
index_i <- floor(seq(1, K, length=seg_i+1))
sortlist1 <- order(rowSums(Data))
sortlist2 <- order(colSums(Data), decreasing=TRUE)

#パラメータとクラス割当の初期値を設定
oldtheta <- matrix(0, nrow=seg_u, ncol=seg_i)
z1 <- rep(0, N)
z2 <- rep(0, K)

for(i in 1:(length(index_u)-1)){
  #ユーザーのクラス割当のインデックスを設定
  index1 <- sortlist1[index_u[i]:index_u[i+1]]
  z1[index1] <- i
  
  for(j in 1:(length(index_i)-1)){
    #アイテムのクラス割当のインデックスを設定
    index2 <- sortlist2[index_i[j]:index_i[j+1]]
    z2[index2] <- j
    
    #セグメントごとのパラメータを決定
    x <- Data[index1, index2]
    oldtheta[i, j] <- sum(x) / length(x)
  }
}

##パラメータの格納用配列
THETA <- array(0, dim=c(seg_u, seg_i, R/keep))
SEG1 <- matrix(0, nrow=R/keep, ncol=N)
SEG2 <- matrix(0, nrow=R/keep, ncol=K)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()

##MCMC推定用配列
LLi1 <- matrix(0, nrow=N, ncol=seg_u)
LLi2 <- matrix(0, nrow=K, ncol=seg_i)
index1 <- list()
index2 <- list()

##データの設定
u_vec <- rep(1, seg_u)
i_vec <- rep(1, seg_i)


####マルコフ連鎖モンテカルロ法でパラメータをサンプリング####
for(rp in 1:R){
  
  ##ユーザーのセグメント割当を生成
  #パラメータを対数変換
  log_theta1 <- log(oldtheta)
  log_theta0 <- log(1-oldtheta)
  
  #セグメントごとの二項分布の対数尤度を計算
  LLi1 <- as.matrix(sparse_data1 %*% t(log_theta1[, z2]) + sparse_data0 %*% t(log_theta0[, z2]))

  #セグメント割当確率のパラメータを設定
  r_matrix <- matrix(r1, nrow=N, ncol=seg_u, byrow=T)   #混合率
  expl <- r_matrix * exp(LLi1 - rowMaxs(LLi1))
  z1_rate <- expl / rowSums(expl)   #セグメント割当確率
  
  #多項分布からセグメント割当を生成
  Zi1 <- rmnom(N, 1, z1_rate)
  z1 <- as.numeric(Zi1 %*% 1:seg_u)
  
  #混合率を更新
  z_sums <- colSums(Zi1) + alpha1
  r1 <- as.numeric(extraDistr::rdirichlet(1, z_sums))
  
  
  ##アイテムのセグメント割当を生成
  #セグメントごとの二項分布の対数尤度を計算
  LLi2 <- as.matrix(sparse_data_T1 %*% log_theta1[z1, ] + sparse_data_T0 %*% log_theta0[z1, ])
  
  #セグメント割当確率のパラメータを設定
  r_matrix <- matrix(r2, nrow=K, ncol=seg_i, byrow=T)   #混合率
  expl <- r_matrix * exp(LLi2 - rowMaxs(as.matrix(LLi2)))
  z2_rate <- expl / rowSums(expl)   #セグメント割当確率
  
  #多項分布からセグメント割当を生成
  Zi2 <- rmnom(K, 1, z2_rate)
  z2 <- as.numeric(Zi2 %*% 1:seg_i)
  
  #混合率を更新
  z_sums <- colSums(Zi2) + alpha2
  r2 <- as.numeric(extraDistr::rdirichlet(1, z_sums))
  
  
  ##ベータ分布からパラメータをサンプリング
  #インデックスを作成
  for(i in 1:seg_u){index1[[i]] <- which(z1==i)}
  for(j in 1:seg_i){index2[[j]] <- which(z2==j)}
  
  #セグメント割当ごとにパラメータを発生させる
  for(i in 1:seg_u){
    for(j in 1:seg_i){
      
      #パラメータ更新に必要なデータを抽出
      x <- Data[index1[[i]], index2[[j]]]   #割当セグメントに該当するデータを抽出
      y <- sum(x)
      n <- length(x)
      
      #ベータ分布からパラメータを発生
      alpha <- tau[1] + y
      beta <- tau[2] + n - y
      oldtheta[i, j] <- rbeta(1, alpha, beta)
    }
  }
  
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta
    SEG1[mkeep, ] <- z1
    SEG2[mkeep, ] <- z2
    
    if(rp%%disp==0){
      LL <- sum((LLi1 * Zi1) %*% u_vec)   #対数尤度
      print(rp)
      print(LL)
      print(round(rbind(r1, mix1), 3))
      print(round(rbind(r2, mix2), 3))
      print(round(cbind(oldtheta, theta0), 3))
    }
  }
}

####サンプリング結果の要約と可視化####
burnin <- 1000/keep   #バーンイン期間 
r <- R/keep   #サンプリングの最終行

##サンプリング結果の可視化
matplot(t(THETA[1, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="二項分布のパラメータ")
matplot(t(THETA[2, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="二項分布のパラメータ")
matplot(t(THETA[3, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="二項分布のパラメータ")
matplot(t(THETA[4, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="二項分布のパラメータ")
matplot(t(THETA[5, , ]), type="l", xlab="サンプリング回数", ylab="パラメータ", main="二項分布のパラメータ")

##パラメータの推定値
round(cbind(apply(THETA[, , burnin:r], c(1, 2), mean), theta0), 3)   #事後平均
round(apply(THETA[, , burnin:r], c(1, 2), function(x) quantile(x, 0.025)), 3)   #事後5％分位点
round(apply(THETA[, , burnin:r], c(1, 2), function(x) quantile(x, 0.975)), 3)   #事後95％分位点
round(apply(THETA[, , burnin:r], c(1, 2), sd), 3)   #事後標準偏差

##セグメント割当と共起関係を可視化
#ユーザーセグメントの割当を確定
seg1 <- rep(0, N)
for(i in 1:N){
  print(i)
  x <- table(SEG1[burnin:r, i])
  seg1[i] <- as.numeric(names(which.max(x)))
}

#アイテムセグメントの割当を確定
seg2 <- rep(0, K)
for(i in 1:K){
  print(i)
  x <- table(SEG2[burnin:r, i])
  seg2[i] <- as.numeric(names(which.max(x)))
}

#データをセグメント順位並び替え
index_seg1 <- order(seg1)
index_seg2 <- order(seg2)
Block_Data <- Data[index_seg1, index_seg2]   #生データをセグメント順に並び替える

#並び替えたブロックデータを可視化
plot_matrix(Data, standardize.rows=FALSE, reorder.rows=FALSE, reorder.cols=FALSE, high.contrast=TRUE)
plot_matrix(Block_Data, standardize.rows=FALSE, reorder.rows=FALSE, reorder.cols=FALSE, high.contrast=TRUE)


