#####階層pitman-yor過程言語モデル#####
options(warn=2)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(data.table)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5723)

####データの発生####
#データの設定
d <- 1000   #文書数
v <- 1000   #語彙数
w <- rpois(d, rgamma(d, 80, 0.40))   #1文書あたりの単語数
f <- sum(w)

#IDの設定
u_id <- rep(1:d, w)
t_id <- c()
for(i in 1:d){
  t_id <- c(t_id, 1:w[i])
}
doc_list <- list()
for(i in 1:d){doc_list[[i]] <- which(u_id==i)}


##bi-gramモデルのパラメータを設定
alpha0 <- sort(rbeta(v, 0.65, 7.0), decreasing=TRUE)
theta01 <- thetat01 <- extraDistr::rdirichlet(1, alpha0)
theta02 <- thetat02 <- extraDistr::rdirichlet(v, alpha0) 
summary(as.numeric(theta02))

##単語を生成
wd0 <- rep(0, f)
WX0 <- matrix(0, nrow=f, ncol=v)
for(i in 1:d){
  if(i%%100==0){
    print(i)
  }
  index <- doc_list[[i]]
  freq <- length(index)
  x <- rep(0, freq)
  
  for(j in 1:freq){
    if(j==1){
      z <- rmnom(1, 1, theta01)
      x[j] <- as.numeric(z %*% 1:v)
      WX0[index[j], ] <- z
      
    } else {
      z <- rmnom(1, 1, theta02[x[j-1], ])
      x[j] <- as.numeric(z %*% 1:v)
      WX0[index[j], ] <- z
    }
  }
  wd0[index] <- x
}
word_freq <- as.numeric((table(c(wd0, 1:v))-1))   #単語の出現数
word_probs <- word_freq / f   #単語の出現確率
hist(word_probs, breaks=25, col="grey", xlab="出現確率", main="単語の出現確率の分布")

#頻度がゼロの単語は除く
index_zeros <- which(colSums(WX0)==0)
WX <- WX0[, -index_zeros]
v <- ncol(WX)
wd <- as.numeric(WX %*% 1:v)
theta1 <- thetat1 <- theta01[, -index_zeros] / sum(theta01[, -index_zeros])
theta2 <- thetat2 <- theta02[-index_zeros, -index_zeros] / rowSums(theta02[-index_zeros, -index_zeros])

##Bi-gramのデータ設計
Data <- matrix(0, nrow=v+1, ncol=v) 
ngram <- matrix(0, nrow=f, ncol=2)
for(i in 1:d){
  index <- doc_list[[i]]
  freq <- length(index)
  for(j in 1:freq){
    if(j==1){
      Data[v+1, wd[index[j]]] <- Data[v+1, wd[index[j]]] + 1
      ngram[index[j], ] <- c(v+1, wd[index[j]])
    } else {
      Data[wd[index[j-1]], wd[index[j]]] <- Data[wd[index[j-1]], wd[index[j]]] + 1
      ngram[index[j], ] <- c(wd[index[j-1]], wd[index[j]])
    }
  }
}

####ギブスサンプリングでBi-gramモデルのパラメータを推定####
##文脈木の初期値を設定
#ユニグラム分布の初期値
n1_data <- rowSums(Data)

#バイグラム分布の初期値
max_table <- 100
index_w2 <- list()
n2_list <- list()
for(i in 1:nrow(Data)){
  index <- which(Data[i, ] > 0)
  n2_list[[i]] <- Data[i, index]
  index_w2[[i]] <- index
}

n2_vec <- unlist(n2_list)
n2_data <- cbind(n2_vec, matrix(0, nrow=length(n2_vec), ncol=max_table-1))
colnames(n2_data) <- 1:max_table
storage.mode(n2_data) <- "integer"


#インデックスを作成
index_n2 <- list()
n <- 0
for(i in 1:nrow(Data)){
  index_n2[[i]] <- n + 1:length(index_w2[[i]])
  n <- max(index_n2[[i]])
}

index_ngram <- ngram
for(i in 1:(v+1)){
  index <- which(index_ngram[, 1]==i)
  for(j in 1:length(index)){
    index_ngram[index[j], 2] <- which(index_w2[[i]]==ngram[index[j], 2])
  }
}


##pitman-yor過程のパラメータを設定
alpha0 <- 0.5
beta0 <- 0.5
G0 <- 1/(v+1)

#
w1_data <- n1_data
w2_data <- n2_data
sum(w2_data)

##単語ごとに文脈木の配置を更新
for(rp in 1:100){
  print(rp)
  for(i in 1:f){
    if(i%%10000==0){
      print(i)
    }
    index1 <- index_n2[[index_ngram[i, 1]]]
    index2 <- index_ngram[i, 2]
    x0 <- w2_data[index1, , drop=FALSE][index2, ]
    if(sum(x0)==1){
      next
    }
    
    ##remove customers
    #テーブルについている客数に比例して客を1人削除
    index <- which(x0 > 0)
    x <- x0[index]
    pr <- x / sum(x)   
    z1 <- as.numeric(rmnom(1, 1, pr))
    x0[index[z1]] <- x[z1]-1
    
    ##pitman-yor過程の基底関数を更新
    #テーブルについている客がいなくなれば親レストランの客も1人削除
    sum(w1_data)
    if(x0[index[z1]]==0){
      w1_data[index2] <- w1_data[index2]-1
      N <- sum(w1_data)
      G1 <- (w1_data-beta0)/(N-1+alpha0) + (alpha0+(v+1)*beta0)/(N-1+alpha0)*G0
    } else {
      N <- sum(w1_data)
      G1 <- (w1_data-beta0)/(N+alpha0) + (alpha0+(v+1)*beta0)/(N+alpha0)*G0
    }
    
    ##add customers
    #テーブルについている客数に比例して客を1人追加
    n <- n2_data[index1, ]
    t1 <- sum(n > 0)
    par1 <- x-length(x)*beta0
    if(min(par1) < 0){
      par1[par1 < 0] <- 0
    }
    par2 <- (alpha0+t1*beta0) * G1[index2]
    pr0 <- par1/(sum(par1)+par2)
    pr <- c(pr0, 1-sum(pr0))
    z2 <- as.numeric(rmnom(1, 1, pr) %*% 1:length(pr))
    x0[z2] <- x0[z2]+1
    
    #新しいテーブルができたら親レストランにも客を1人追加
    if(x0[z2]==1){
      w1_data[index2] <- w1_data[index2]+1
    }
    
    if(length(index1)==1){
      w2_data[index1, ] <- x0
    } else {
      w2_data[index1, ][index2, 1:length(x0)] <- x0
    }
  }
  print(sum(w2_data))
}
sum(w1_data)
w2_data
cbind(rowSums(w2_data), rowSums(w2_data > 0))

w2_data[, 1:20]

cbind(G1, w1_data/sum(w1_data))


sum(w2_data)
sum(n2_data)
