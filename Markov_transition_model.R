#####マルコフモデルによるWebサイト回遊行動分析#####
library(MASS)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
#set.seed(93417)
N <- 2000   #サンプル数
S <- 5   #サイトページ数

##サイト回遊確率と離脱確率の設定
#マルコフ推移行列の設定
p1 <- c(0.45, 0.2, 0.1, 0.05, 0.2)
p2 <- c(0.3, 0.3, 0.1, 0.2, 0.1)
p3 <- c(0.3, 0.2, 0.2, 0.1, 0.2)
p4 <- c(0.3, 0.3, 0.1, 0.15, 0.15)
p5 <- c(0.2, 0.3, 0.1, 0.2, 0.2)
p <- rbind(p1, p2, p3, p4, p5)   #データの結合
pt <- p

#ランディングページの確率の設定
pf <- c(0.5, 0.3, 0.1, 0.05, 0.05)
pft <- pf

#離脱確率の関数の設定
surv <- function(beta, X){
  logit <- beta[1] + as.matrix(X) %*% beta[-1]
  pr <- exp(logit)/(1+exp(logit))
  return(pr)
}

#離脱モデルの回帰係数の設定
alpha <- -2.0
beta.time <- 0.10
beta.page <- c(-0.6, -0.4, 0.15, 0.30)
beta <- c(alpha, beta.time, beta.page)
betat <- beta

x <- t(rmultinom(30, 1, runif(5)))
X <- cbind(1:30, x[, -3])
surv(beta, X)


##離散時間マルコフ生存モデルでシミュレーションデータを発生
#データの格納用オブジェクトの設定
X <- matrix(0, 200000, S+1)
Z <- c()
PZ <- c()

for(i in 1:N){
#ランディングページを決定
  print(i)
  r <- sum(rowSums(X)!=0)
  l <- t(rmultinom(1, 1, pf))
  x <- t(c(1, l[-3])) 
  X[r+1, ] <- t(c(1, l))
  
  #2回目以降サイト回遊
  M <- matrix(0, 100, S)
  
  for(t in 2:1000){
    #離脱したかどうかを決定
    pz <- surv(beta=beta, X=x)   #離脱確率
    z <- rbinom(1, 1, pz)   #離脱したかどうかの決定
    PZ <- c(PZ, pz)
    Z <- c(Z, z)
    if(z==1) break   #z=1で離脱
    
    #どのページに進んだかを決定
    l <- (1-z) * (rmultinom(1, 1, p[which.max(l), ]))   #マルコフ推移
    M[1, ] <- t(l)
    x <- t(c(t, l[-3]))
    X[r+t, ] <- t(c(t, l))
  }
}

##不要なデータを削除する
index <- rowMeans(X)!=0   #余っているデータ部分を特定
X <- X[index, ]
summary(X)   #データの要約

##IDの設定
id <- c()
for(i in 1:N){
  r <- rep(i, subset(1:length(Z), Z==1)[i]-length(id))
  id <- c(id, r)  
}

##すべてのデータを結合
ZX <- data.frame(id, Pr=PZ, Z, t=X[, 1], X=X[, -1])
round(ZX, 2)


####最尤法で離散時間生存マルコフ推移モデルを推定####
#マルコフ性が成り立つので、ランディングページ、ページ遷移、離脱はそれぞれ独立と仮定できる。
#したがって、3つのパラメータは独立に求める。

##ページそれぞれのランディング確率を求める
Pl <- colMeans(X[X[, 1]==1, 2:ncol(X)])
round(Pl, 3)


##ページ推移確率行列を求める
#ランディングページで離脱しているユーザーを取り除く
index1 <- subset(1:N, as.numeric(table(id))!=1)
R <- c()
for(i in 1:length(index1)){
  r <- subset(1:length(id), id==index1[i])
  R <- c(R, r)
}
Xm <- X[R, ]   #ページ推移があるデータのみ取り出す

#マルコフ推移行列を求める
Pr <- matrix(0, 0, S) 
for(s in 1:S){
  index1 <- subset(1:nrow(Xm), Z[R]!=1 & Xm[, s+1]==1)   #Z=1のデータを取り除き、t-1のデータを取り出す
  Xm1 <- Xm[index1, 2:ncol(Xm)]   #tの列を取り除く
  Nm <- sum(Xm1[, s])   #t-1でsのデータを合計する
  Mm <- colSums(Xm[index1+1, 2:ncol(Xm)])   #tのデータを合計する
  pr <- Mm / Nm   #確率を計算
  Pr <- rbind(Pr, pr)
}
rownames(Pr) <- c("P1", "P2", "P3", "P4", "P5")
round(Pr, 3)


##離散時間生存モデルで離脱率を推定
#離散時間生存モデルの対数尤度を設定
loglike <- function(beta, y, X){
  logit <- beta[1] + as.matrix(X) %*% beta[-1]    #ロジットの計算
  p <- exp(logit)/(1+exp(logit))   #確率の計算
  
  #対数尤度の計算
  LLs <- y*log(p) + (1-y)*log(1-p)
  LL <- sum(LLs)
  return(LL)
}

#対数尤度を最大化する
b0 <- c(-1.0, 0.2, runif(S-1, -1, 1))
res <- optim(b0, loglike, y=Z, X=X[, -4], method="BFGS", hessian=TRUE, control=list(fnscale=-1))

#推定結果
beta <- res$par
round(beta, 3)   #推定された回帰係数
round(betat, 3)   #真の回帰係数


##すべての推定結果を表示
#左が推定結果、右が真の値
round(Pl, 2); round(pft, 3)   #ページ別のランディング確率
round(Pr, 2); round(pt, 2)   #マルコフ推移行列
round(beta, 2); round(betat, 2)   #サイト離脱モデルのパラメータ


####離脱率のシミュレーション####
n <- 500   #シミュレーションするサンプル数
t <- 100   #シミュレーションする期間

#データの格納用オブジェクトの設定
Xs <- matrix(0, n*t, S+1)
Zs <- c()
PZs <- c()
for(i in 1:n){
  #ランディングページを決定
  print(i)
  r <- sum(rowSums(Xs)!=0)
  l <- t(rmultinom(1, 1, Pl))
  x <- t(c(1, l[-3])) 
  Xs[r+1, ] <- t(c(1, l))
  
  #2回目以降サイト回遊
  M <- matrix(0, 100, S)
  
  for(t in 2:t){
    #離脱したかどうかを決定
    p <- surv(beta=beta, X=x)   #離脱確率
    z <- rbinom(1, 1, p)   #離脱したかどうかの決定
    PZs <- c(PZs, p)
    Zs <- c(Z, z)
    
    #どのページに進んだかを決定
    l <- (rmultinom(1, 1, Pr[which.max(l), ]))   #マルコフ推移
    M[1, ] <- t(l)
    x <- t(c(t, l[-3]))
    Xs[r+t, ] <- t(c(t, l))
  }
  #最終時間の離脱率
  p <- surv(beta=beta, X=x)   #離脱確率
  z <- rbinom(1, 1, p)   #離脱したかどうかの決定
  PZs <- c(PZs, p)
  Zs <- c(Z, z)
}

#ページ推移数のみを考慮した離脱率
s <- matrix(0, t, S-1)
Xt <- data.frame(1:t, s)
X.surv <- surv(beta=beta, X=Xt)

##データを結合
ids <- rep(1:n, rep(t, n))   #idを設定
ZXs <- data.frame(id=ids, p=PZs, t=Xs[, 1], S=Xs[, -1])   #データを結合
round(ZXs, 3)


##離脱率を可視化
plot(1:t, X.surv, type="l", ylab="離脱率", xlab="ページ推移数", lwd=3, col=2, ylim=c(0.1, 1))
for(i in 1:5){
  lines(1:t, ZXs[ZXs$id==i, ]$p, type="l", lty=i)
}
lines(1:t, X.surv, type="l", ylab="離脱率", xlab="ページ推移数", lwd=3, col=2, ylim=c(0.1, 1))
