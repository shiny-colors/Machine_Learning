#####k-means法######
library(MASS)
library(dplyr)
library(reshape2)
####データの発生#####
#set.seed(493)
k <- 4   #クラスター数
val <- 6   #説明変数の数
n <- 500   #セグメントごとのデータ数

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
corM <- corrM(col=val, lower=-0.2, upper=0.2)
eigen(corM)

##相関行列から分散共分散行列を作成する関数を定義
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
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

#分散共分散行列を作成
Sigma1 <- covmatrix(col=val, corM=corM, lower=25, upper=35)
Sigma2 <- covmatrix(col=val, corM=corM, lower=30, upper=40)
Sigma3 <- covmatrix(col=val, corM=corM, lower=20, upper=30)
Sigma4 <- covmatrix(col=val, corM=corM, lower=27, upper=38)
Sigma <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance, Sigma3$covariance)

#群ごとの変数の平均を作成
mu1 <- c(rnorm(val, 14, 10))
mu2 <- c(rnorm(val, 19, 10))
mu3 <- c(rnorm(val, 9, 8))
mu4 <- c(rnorm(val, 23, 13))
mu <- list(mu1, mu2, mu3, mu4)

##多変量正規分布からの乱数を発生させる
k; n
X <- matrix(0, 0, val)
for(kk in 1:k){
  xx <- mvrnorm(n=n, mu[[kk]], Sigma[[kk]])
  X <- rbind(X, xx)
}

#クラスターのベクトルを作成
y <- rep(1:4, rep(n, 4))
table(y)

#データを結合してデータを要約
YX <- as.data.frame(cbind(y, X))
by(YX, YX[, 1] , function(x) summary(x))   #セグメントごとの要約関数
by(YX[, 2:7], YX[, 1], colMeans)   #セグメントごとの平均
by(YX[, 2:7], YX[, 1] , var)   #セグメントごとの分散共分散行列

####k-means法でセグメントに分ける####
##初期分割を設定
#変数を標準化して行の合計値を昇順に並べて初期分割とする
Xscale <- scale(YX[, -1])
Xsums <- cbind(YX, apply(Xscale, 1, sum))
names(Xsums)[8] <- "Xscale"
Xseg <- Xsums[order(Xsums[, 8], decreasing = TRUE), c(-1, -8)]
trueSeg <- Xsums[order(Xsums[, 8], decreasing = TRUE), 1]   #真のセグメントも順番に並べておく
ns <- dim(Xseg)[1]/k   #1セグメントあたりの分割数

#セグメントごとに分割して平均ベクトルを求める
X1mean <- colMeans(Xseg[1:ns, ])
X2mean <- colMeans(Xseg[(ns+1):(ns*2), ])
X3mean <- colMeans(Xseg[(ns*2+1):(ns*3), ])
X4mean <- colMeans(Xseg[(ns+1):(ns*4), ]) 
Xmean <- list(X1mean, X2mean, X3mean, X4mean)

##それぞれのサンプルに対してセグメントごとのユークリッド距離を求めて、
##その距離が最小となるセグメントにサンプルを所属させ初期分割を作る
segsign <- rep(1:4, rep(500, 4))
Sw <- numeric()
for(m in 1:k){
  ss <- apply(Xseg[segsign==m, ], 1, function(x) sum((x-Xmean[[m]])^2))
  (ssm <- sum(ss)/table(segsign)[[m]])
  Sw <- c(Sw, ssm)  
}
(Sold <- sum(Sw))   #クラスター平方和の合計を計算
diff <- 100   #平方和の差の初期値
tol <- 1   #停止基準


###k-means法のアルゴリズム
while(abs(diff) >= tol){
  ##所属セグメントを決定するアルゴリズム
  Xdis <- matrix(0, nrow=nrow(Xseg), ncol=k)
  for(kk in 1:k){
    xxd <- apply(Xseg, 1, function(x) sum((x-Xmean[[kk]])^2))
    Xdis[, kk] <- xxd
  }
  #所属するセグメントの設定
  seg <- apply(Xdis, 1, which.min)   #平方和が最小のセグメントに所属させる
  segN <- table(seg)
  YXseg <- cbind(seg, Xseg)
  
  ##セグメント分割の評価を求めるアルゴリズム
  #セグメントごとに平均ベクトルを求める
  Xmean <- list()
  for(m in 1:k){
    xm <- apply(Xseg[seg==m, ], 2, mean)
    Xmean[[m]] <- xm
  }
  
  #クラスター内平方和を計算
  Sw <- numeric()
  for(m in 1:k){
    ss <- apply(Xseg[seg==m, ], 1, function(x) sum((x-Xmean[[m]])^2))
    (ssm <- sum(ss)/segN[m])
    Sw <- c(Sw, ssm)  
  }
  S <- sum(Sw)   #クラスター平方和の合計を計算
  diff <- abs(Sold - S)   #以前のクラスター平方和との差の絶対値を計算
  Sold <- S
  print(S)
}

##結果の要約
resultX <- cbind(trueSeg, seg, Xseg)   #真のセグメント、推定されたセグメント、データを結合
sortlist <- order(resultX[, 2])
resultX <- resultX[sortlist, ]
rownames(resultX) <- c(1:nrow(resultX))   #行番号を振り直す 
table(resultX[, 2], resultX[, 1])   #真のセグメントと推定されたセグメントの誤判別表

#推定されたセグメントごとの要約
by(resultX[, 3:8], resultX[, 2] , function(x) summary(x))   #推定されたセグメントごとの要約関数
by(resultX[, 3:8], resultX[, 2], colMeans)   #セグメントごとの平均
by(resultX[, 3:8], resultX[, 2] , var)   #セグメントごとの分散共分散行列
plot(resultX[, 3:8], col=resultX[, 2])

#真のセグメントごとの要約
by(resultX[, 3:8], resultX[, 1] , function(x) summary(x))   #推定されたセグメントごとの要約関数
by(resultX[, 3:8], resultX[, 1], colMeans)   #セグメントごとの平均
by(resultX[, 3:8], resultX[, 1] , var)   #セグメントごとの分散共分散行列
plot(resultX[, 3:8], col=resultX[, 1])

####関数を使う####
res <- kmeans(x=resultX[, 3:8], centers=k)
resX <- cbind(res$cluster, resultX)
names(resX)[1] <- "cluster_f"
table(resX$seg, resX$cluster_f)   #自作のkmeans法とRのkmeans法を比較


