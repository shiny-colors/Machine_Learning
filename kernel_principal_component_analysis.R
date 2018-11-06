#####主成分分析#####
library(MASS)
library(plyr)
library(reshape2)
library(kernlab)

####データの発生####
#set.seed(4238)
#多変量正規分布からの乱数発生
val <- 10   #説明変数の数
n <- 1000   #サンプル数

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
corM <- corrM(col=10, lower=0.3, upper=0.9)
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
Sigma <- covmatrix(col=10, corM=corM, lower=20, upper=30)

#群ごとの変数の平均を作成
mu <- c(rnorm(10, 20, 15))

##多変量正規分布からの乱数を発生させる
val; n
X <- mvrnorm(n=n, mu, Sigma$covariance)

##データを要約
round(colMeans(X), 2)   #変数ごとの平均
round(var(X), 2)   #分散共分散行列
round(cor(X), 2)   #相関行列
summary(X)   #要約

##散布図行列の作成
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#変数1～5の散布図行列
pairs(as.data.frame(X[, 1:5]), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
      upper.panel=panel.cor)
#変数6～10の散布図行列
pairs(as.data.frame(X[, 6:10]), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
      upper.panel=panel.cor)

####主成分分析を実行####
##生データから主成分分析
Xvar <- var(X)   #分散共分散行列を求める
PCA <- eigen(Xvar)   #固有値分解(主成分分析)を行う
PCAval <- PCA$values   #固有値
PCAvec <- PCA$vectors   #固有ベクトル
round(PCAvec %*% diag(PCAval) %*% t(as.matrix(PCAvec)) - Xvar, 3)   #特異値分解が出来ているか見る

#PCAの要約
round(PCAval/sum(PCAval), 3)   #成分ごとの寄与率
round(cumsum(PCAval)/sum(PCAval), 3)   #累積寄与率

#第3主成分まで主成分得点を計算
PC1 <- X %*% PCAvec[, 1]
PC2 <- X %*% PCAvec[, 2]
PC3 <- X %*% PCAvec[, 3]
PC <- data.frame(PC1, PC2, PC3)
plot(PC, col=4)

####カーネル主成分分析####
##グラム行列の作成
#多項式カーネル
Xscale <- scale(X)
gram1 <- (1+as.matrix(X) %*% t(as.matrix(X))) + (1+as.matrix(X) %*% t(as.matrix(X)))^2
round(gram1[1:15, 1:15], 3)
round(gram1[985:1000, 985:1000], 3)
round(eigen(gram1)$value, 3)   #半正定値がどうか確認

#ガウスカーネル
sigma <- 1/2
kf_gauss <- function(x1, sigma){
  x1 <- as.matrix(x1)
  g1 <- matrix(t(x1), nrow(x1)^2, ncol(x1), byrow=T)
  g2 <- matrix(rep(x1, c(rep(nrow(x1), nrow(x1)*ncol(x1)))), nrow(x1)^2, ncol(x1))
  g3 <- (g2 - g1)^2
  gram <- exp(-sigma*sqrt(matrix(apply(g3, 1, sum), nrow(x1), nrow(x1), byrow=T)))
  return(gram)
}
gram2 <- kf_gauss(x1=X, sigma=sigma)
round(gram2[1:5, 1:5], 3)
round(gram2[985:1000, 985:1000], 3)
round(eigen(gram2)$value, 3)   #半正定値がどうか確認

#関数を使う
L <- kernelMatrix(rbfdot(sigma=1/2), X)   #ガウスカーネルで変換

##グラム行列から主成分分析を実行
PCA <- eigen(gram1)   #固有値分解(主成分分析)を行う
PCAval <- PCA$values   #固有値
PCAvec <- PCA$vectors   #固有ベクトル
round(PCAvec %*% diag(PCAval) %*% t(as.matrix(PCAvec)) - gram1, 3)   #特異値分解が出来ているか見る

#PCAの要約
round(PCAval/sum(PCAval), 3)   #成分ごとの寄与率
round(cumsum(PCAval)/sum(PCAval), 3)   #累積寄与率

#第3主成分まで主成分得点を計算
PC1 <- gram1 %*% PCAvec[, 1]
PC2 <- gram1 %*% PCAvec[, 2]
PC3 <- gram1 %*% PCAvec[, 3]
PC <- data.frame(PC1, PC2, PC3)
plot(PC, col=4)
 
##関数を使う
kp <- kpca(X, kernel="polydot", kpar=list(degree=2, scale=1))
plot(data.frame(pcv(kp))[, 1:3])
