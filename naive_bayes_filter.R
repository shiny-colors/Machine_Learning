#####ナイーブベイズフィルタ#####
library(MASS)
library(e1071)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
Nnm <- 1500   #正常データのサンプル数
Nab <- 1500   #異常データのサンプル数
N <- 3000   #総サンプル数
k <- 200   #語彙数
w <- 500   #文書あたりの単語数
z <- c(rep(0, Nnm), rep(1, Nab))   #正常と異常のインディケーター

#分類に関連しないデータの分布の発生
W1 <- matrix(0, N, 3/4*k) 
for(i in 1:(3/4*k)){
  p1 <- runif(1, 0.01, 0.25)
  W1[, i] <- rbinom(N, 1, p1)
}

#正常データに関連するデータの発生
Wnm <- matrix(0, N, 1/4*k/2)
for(i in 1:(1/4*k/2)){
  pnm1 <- runif(1, 0.075, 0.3)
  pnm2 <- runif(1, 0.005, 0.1)
  Wnm1 <- rbinom((N/2), 1, pnm1)
  Wnm2 <- rbinom((N/2), 1, pnm2)
  Wnm.r <- c(Wnm1, Wnm2)
  Wnm[, i] <- Wnm.r
}

#異常データに関連するデータの発生
Wab <- matrix(0, N, 1/4*k/2)
for(i in 1:(1/4*k/2)){
  pab1 <- runif(1, 0.01, 0.1)
  pab2 <- runif(1, 0.15, 0.3)
  Wab1 <- rbinom((N/2), 1, pab1)
  Wab2 <- rbinom((N/2), 1, pab2)
  Wab.r <- c(Wab1, Wab2)
  Wab[, i] <- Wab.r
}

#データの結合
Wz <- cbind(z, W1, Wnm, Wab)

#データの集計
by(Wnm, z, colSums)
by(Wab, z, colSums)


####ナイーブベイズ分類器####
#学習データとテストデータに分ける
Wz.learn <- rbind(Wz[1:1000, ], Wz[1501:2500, ])
Wz.test <- rbind(Wz[1001:1500, ], Wz[2501:3000, ])
z.learn <- c(z[1:1000], z[1501:2500])
z.test <- c(z[1001:1500], z[2501:3000])

##ナイーブベイズの尤度関数
#単語頻度を数える
wnm <- colSums(Wz.learn[z.learn==0, -1])   #正常データでの単語頻度
wab <- colSums(Wz.learn[z.learn==1, -1])   #異常データでの単語頻度

#条件付き確率
round(pnm <- wnm / (wnm + wab), 3)   #正常データの条件付き確率
round(pab <- wab / (wnm + wab), 3)   #異常データの条件付き確率

#事前確率　
prior <- 0.5

##ベイズフィルタで正常データ、異常データを判定(事前確率は同じ　)
score <- c()
for(i in 1:nrow(Wz.test)){
  #正常、異常の尤度
  post.nm <- Wz.test[i, -1]*pnm
  post.ab <- Wz.test[i, -1]*pab
  
  #尤度がゼロ以外のデータを取り出す
  post.nmz <- subset(post.nm, post.nm > 0)
  post.abz <- subset(post.ab, post.ab > 0)
  
  #スコアを計算
  score <- c(score, sum(log(post.abz))-sum(log(post.nmz)))
}

#結果
round(score, 3)
(res <- ifelse(score > 0, 1, 0))
(res.table <- table(z.test, res))   #誤判別表
c(res.table[1, 1] / sum(res.table[1, ]), res.table[2, 2] / sum(res.table[2, ]))   #分類別の正答率
sum(diag(res.table)) / sum(res.table)   #全体での正答率


