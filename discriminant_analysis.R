#####判別モデル#####
library(MASS)
library(plyr)
library(reshape2)
####多群の正準判別モデル####
####データの発生####
#set.seed(4235)
#多変量正規分布からの乱数発生
k <- 4   #群の数
val <- 6   #説明変数の数
n <- 400   #群ごとの学習に使うためのデータ数
nt <- 300   #群ごとのテストに使うためのデータ数

##多変量正規分布からの乱数を発生させる
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#群ごとの相関行列を作成(群ですべて同じ)
corM <- corrM(col=6, lower=-0.2, upper=0.2)
eigen(corM)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
    diag(c) <- m   #対角行列を元の分散に戻す
  }
  diag(c) <- m
  cc <- c * corM
  #固有値分解で強制的に正定値行列に修正する
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

#分散共分散行列を作成(群ですべて同じ)
Sigma1 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma2 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma3 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma4 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
S <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance, Sigma4$covariance)

#群ごとの変数の平均を作成
mu1 <- c(rnorm(6, 12, 10))
mu2 <- c(rnorm(6, 18, 10))
mu3 <- c(rnorm(6, 6, 10))
mu4 <- c(rnorm(6, 24, 10))
mu <- list(mu1, mu2, mu3, mu4)

##多変量正規分布からの乱数を発生させる
k; n; nt
X <- matrix(0, 0, 6)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], S[[kk]])
  X <- rbind(X, xx)
}

#教師データのベクトルを作成
y <- rep(1:4, rep(700, 4))

#データを結合
YX <- as.data.frame(cbind(y, X))
by(YX, YX[, 1] , function(x) summary(x))   #群ごとの要約関数
by(YX[, 2:7], YX[, 1] , function(x) cor(x))   #群ごとの相関
plot(YX[, 2:7], col=YX[, 1])   #散布図

#テストデータと学習データを分離
#学習データ
YXl <- rbind(YX[1:400, ], YX[701:1100, ], YX[1401:1800, ], YX[2101:2500, ])
table(YXl[, 1])

#テストデータ
YXt <- rbind(YX[401:700, ], YX[1101:1400, ], YX[1801:2100, ], YX[2501:2800, ])
table(YXt[, 1])

####正準判別モデルで学習データから分類器を学習する####
#群内平均
mucat <- matrix(0, k, val)
for(i in 1:k){
  muc <- apply(YXl[YXl[, 1]==i, 2:7], 2, mean)
  mucat[i, ] <- (muc)
}

#群全体の平均
muall <- colMeans(YXl[, 2:7])
names(muall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")

#群全体の分散共分散行列
covall <- var(YXl[, 2:7])
rownames(covall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")
colnames(covall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")

#群間の分散共分散行列
covB <- 1/k * t(mucat - muall) %*% (mucat - muall)

##固有値問題を解いて群間の分離度を最大化する解を得る
covA <- solve(covall)
M <- eigen(covA %*% as.matrix(covB))   #固有値問題を解く
R <- M$values   #行列のランクは群数-1
a <- M$vectors   #固有ベクトルが判別関数の係数
a
(cont <- R/sum(R))   #寄与率
(cumcont <-cumsum(cont))   #累積寄与率

#第2固有値までを用いて、2次元空間上に射影した行列をプロット
#学習データに対する合成変量とプロット
Zl <- as.matrix(YXl[, 2:7]) %*% t(a[1:3, ])
plot(as.data.frame(Zl), col=YXl[, 1])

#テストデータに対する合成変量とプロット
Zt <- as.matrix(YXt[, 2:7]) %*% t(a[1:3, ])
plot(as.data.frame(Zt))

####判別性能を確認####
#学習データの正答率を計算
zm <- matrix(0, k, 3)
for(i in 1:k){
  zm[i, ] <- apply(Zl[YXl[, 1]==i, ], 2, mean)
}

#判別結果を求める
D <- matrix(0, nrow=nrow(ZZl), ncol=k)
for(i in 1:k){
  ZZl <- Zl - matrix(zm[i, ], nrow=nrow(Zl), ncol=ncol(Zl), byrow=T)
  d <- matrix(apply(ZZl, 1, function(x) t(x) %*% x), nrow=nrow(ZZl), ncol=1)
  D[, i] <- d
}
res <- apply(D, 1, which.min)   #判別結果
pre <- data.frame(correct=YXl[, 1], predict=res)   #正解データと結合
table(pre)   #誤判別表
round(apply(table(pre), 1, function(x) x/sum(x)), 3)   #群ごとの正答率
sum(diag(table(pre)))/sum(table(pre))   #正答率


####2次判別####
##異なる共分散行列の多変量正規分布から乱数発生
#群ごとの相関行列を作成(群で異なる)
k <- 3   #群数
val <- 6   #説明変数数
n <- 500   #群ごとの学習データ
nt <- 500   #群ごとのテストデータ
N <- 3000   #全サンプル数

corM1 <- corrM(col=6, lower=-0.25, upper=0.3)
corM2 <- corrM(col=6, lower=-0.3, upper=0.4)
corM3 <- corrM(col=6, lower=-0.4, upper=0.55)

#分散共分散行列を作成(群ですべて同じ)
Sigma1 <- covmatrix(col=6, corM=corM1, lower=15, upper=26)
Sigma2 <- covmatrix(col=6, corM=corM2, lower=18, upper=30)
Sigma3 <- covmatrix(col=6, corM=corM3, lower=12, upper=18)
S <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance)

#群ごとの変数の平均を作成
mu1 <- c(rnorm(6, 12, 10))
mu2 <- c(rnorm(6, 18, 10))
mu3 <- c(rnorm(6, 8, 10))
mu <- list(mu1, mu2, mu3)

##多変量正規分布からの乱数を発生させる
X <- matrix(0, 0, 6)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], S[[kk]])
  X <- rbind(X, xx)
}

#教師データのベクトルを作成
y <- rep(1:3, rep(1000, 3))

#データを結合
YX <- data.frame(y, X)
by(YX, YX[, 1] , function(x) summary(x))   #群ごとの要約関数
by(YX[, 2:7], YX[, 1] , function(x) cor(x))   #群ごとの相関
plot(YX[, 2:7], col=YX[, 1])   #散布図

#テストデータと学習データを分離
#学習データ
YXl <- rbind(YX[1:500, ], YX[1001:1500, ], YX[2001:2500, ])
table(YXl[, 1])

#テストデータ
YXt <- rbind(YX[501:1000, ], YX[1501:2000, ], YX[2501:3000, ])
table(YXt[, 1])

####2次判別の推定####
##群ごとの平均ベクトルと分散共分散行列を計算
#群ごとの平均ベクトル
mv <- matrix(0, k, val)
for(mm in 1:k){
  m <- apply(YXl[YXl[, 1]==mm, 2:7], 2, mean)
  mv[mm, ] <- m
}

#群ごとの分散共分散行列
Sk <- list()
for(ss in 1:k){
  m <- var(YXl[YXl[, 1]==ss, 2:7])
  Sk[[ss]] <- m
}

##サンプルごとにそれぞれの群のマハラノビスの汎距離を求めてもっとも小さい群に所属させる
##学習データについて求める
HL <- matrix(0, nrow=nrow(YXl), ncol=k)
S <- list(solve(Sk[[1]]), solve(Sk[[2]]), solve(Sk[[3]]))
for(kk in 1:k){
  h <- apply(YXl[, 2:7], 1, function(x) t(x-mv[kk, ]) %*% S[[kk]] %*% (x-mv[kk, ]))
  HL[, kk] <- h
}

#学習データの判別結果
pred <- apply(HL, 1, which.min)   #判別結果
res <- data.frame(true=YXl[, 1], pred)   #真の群と結合
table(res)   #誤判別表
sum(diag(table(res)))/sum(table(res))   #誤判別率


##テストデータについて求める
HT <- matrix(0, nrow=nrow(YXt), ncol=k)
S <- list(solve(Sk[[1]]), solve(Sk[[2]]), solve(Sk[[3]]))
for(kk in 1:k){
  h <- apply(YXt[, 2:7], 1, function(x) t(x-mv[kk, ]) %*% S[[kk]] %*% (x-mv[kk, ]))
  HT[, kk] <- h
}

#テストデータの判別結果
pred <- apply(HT, 1, which.min)   #判別結果
res <- data.frame(true=YXt[, 1], pred)   #真の群と結合
table(res)   #誤判別表
sum(diag(table(res)))/sum(table(res))   #誤判別率