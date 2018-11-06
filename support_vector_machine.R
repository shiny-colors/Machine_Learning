#####サポートベクターマシン#####
library(MASS)
library(kernlab)
library(mlbench)
library(dplyr)
library(reshape2)
library(quadprog)
####データの発生#####
k <- 2   #群の数
val <-12   #説明変数の数
n <- 500   #群ごとの学習に使うためのデータ数
nt <- 500   #群ごとのテストに使うためのデータ数

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

#分散共分散行列を作成(群ですべて同じ)
Sigma1 <- covmatrix(col=val, corM=corM, lower=25, upper=35)
Sigma2 <- covmatrix(col=val, corM=corM, lower=30, upper=40)
Sigma <- list(Sigma1$covariance, Sigma2$covariance)
Sigma
#群ごとの変数の平均を作成
mu1 <- c(rnorm(val, 10, 10))
mu2 <- c(rnorm(val, 14, 10))
mu <- list(mu1, mu2)

##多変量正規分布からの乱数を発生させる
k; n; nt
X <- matrix(0, 0, val)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], Sigma[[kk]])
  X <- rbind(X, xx)
}

#教師データのベクトルを作成
y <- rep(c(1, 0), rep(n+nt, 2))
table(y)

#一部の変数を離散変数に変換
mu_rbind <- colMeans(rbind(mu[[1]], mu[[2]]))   #2群の統合平均
x2 <- t(apply(X[, 8:val], 1, function(x) ifelse(x > mu_rbind[8:val], 1, 0)))   #統合平均より高い要素を1とする
mu_rbind[8:val]
head(cbind(X[, 8:val], x2), 20)
head(YX)

#連続変数を標準化
Xscale <- scale(X[, 1:7])

#データを結合
YX <- as.data.frame(cbind(y, Xscale, x2))
by(YX, YX[, 1] , function(x) summary(x))   #群ごとの要約関数
by(YX[, 2:8], YX[, 1] , function(x) cor(x))   #群ごとの連続変数の相関
by(YX[, 9:val+1], YX[, 1], colMeans)   #群ごとの離散変数の出現率

#連続変数の説明変数をプロット
yplotcolor <- YX[, 1]
yplotcolor[yplotcolor==0] <- 2
plot(YX[, 2:8], col=yplotcolor)   #連続変数の散布図
YX[YX[, 1]==0, 1] <- -1   #0の群を-1にする

#学習データとテストデータに分ける
YXl <- YX[c(1:500, 1001:1500), ]
YXt <- YX[c(501:1000, 1501:2000), ]

####サポートベクターマシンの実装####
##線形サポートベクターマシン
#最大化する関数
dim(YXl)
yy <- matrix(YXl[, 1], nrow=nrow(YXl), ncol=ncol(YXl[, -1])) 
xx <- yy * YXl[, -1]
Q <- as.matrix(xx[, -1]) %*% t(xx[, -1])
c <- matrix(rep(-1, nrow(YXl)))

#制約条件
sm <- 0.0001   #ソフトマージン
A <- t(YXl[, 1])
b <- 0
l <- matrix(rep(0, nrow(YXl)))
u <- matrix(rep(sm), nrow(YXl))
r <- 0

#凸二次計画問題を解く
sv <- ipop(c, Q, A, b, l, u, r)
sv
dual(sv)

#betaを求める
a <- primal(sv)
aa <- a>0   #a>0のみ取り出す
alpha <- matrix(a[aa], nrow=nrow(xx[aa, ]), ncol=ncol(xx))
(beta <- colSums(alpha * xx[aa, ]))   #サポートベクターを用いてbetaを算出

#betaからbを求める
b <- mean(YXl[aa, 1] - as.matrix(YXl[aa, -1]) %*% as.matrix(beta))

##分類結果
#学習データの場合
resl <- as.matrix(YXl[, -1]) %*% as.matrix(beta) + b    #決定境界
resll <- sign(resl)
sum(as.numeric(resll[1:500]==1))/500
sum(as.numeric(resll[501:1000]==-1))/500

#テストデータの場合
rest <- as.matrix(YXt[, -1]) %*% as.matrix(beta) + b   #決定境界
restt <- sign(rest)
sum(as.numeric(restt[1:500]==1))/500
sum(as.numeric(restt[501:1000]==-1))/500

#判別結果のクロス集計
testerror1 <- cbind(YXt[, 1], restt)
table(testerror1[, 1], testerror1[, 2])
 
##関数を使う場合
model <- ksvm(factor(YXl[, 1]) ~ ., data=YXl, kernel="vanilladot", prob.model=FALSE)
model
print(model@alpha)   #alphaの係数
print(model@b)   #betaの係数
print(model@SVindex)   #サポートベクタ-

#判別結果
pre <- predict(model, as.data.frame(YXt), type="response")   #予測結果
pre <- as.numeric(pre)
pre[pre==1] <- -1 
pre[pre==2] <- 1
testerror2 <- cbind(YXt[, 1], pre)
table(testerror2[, 1], testerror2[, 2])

####カーネル法を用いたサポートベクターマシン####
##カーネル関数の定義
##グラム行列を作成する
#多項式カーネル
gram1 <- (3 + as.matrix(Xnew) %*% t(as.matrix(Xnew))) + (3 + as.matrix(Xnew) %*% t(as.matrix(Xnew)))^2
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
gram2 <- kf_gauss(x1=Xnew, sigma=sigma)
round(gram2[1:15, 1:15], 3)
round(gram2[985:1000, 985:1000], 3)
round(eigen(gram2)$value, 3)   #半正定値がどうか確認

#関数を使う
L <- kernelMatrix(rbfdot(sigma=1/2), Xnew)   #ガウスカーネルで変換

