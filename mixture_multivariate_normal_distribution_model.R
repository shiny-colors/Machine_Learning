#####混合正規分布(Smart Shifter)#####
library(MASS)
library(mclust)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####多変量正規分布の乱数を発生させる関数を定義####
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

####データの発生####
#set.seed(421830)
##データの設定
seg <- trunc(runif(5, 2, 6.3))   #セグメント数
N <- 5000   #サンプル数
col <- 8   #説明変数数
ca <- 5   #カテゴリー数

##カテゴリカル変数の発生
p <- runif(ca, 0.2, 1.0) 
category <- apply(t(rmultinom(N, 1, p)), 1, which.max)
sortlist <- order(category)
category <- category[sortlist]
table(category)

##多変量正規乱数を用いて変数を発生させる
#変数とパラメータの格納
Xlist <- list()
MU <- list()
covM <- list()
Z <- c()

##カテゴリーおよびセグメントごとに変数を発生
for(c in 1:ca){
  X.seg <- list()   #リストを初期化

  ##カテゴリごとに混合正規分布を発生
  ps <- runif(seg[c], 0.2, 1.0)
  n.seg <- length(category[category==c])
  seg.z <- apply(t(rmultinom(n.seg, 1, ps)), 1, which.max)
  sort.seg <- order(seg.z)
  seg.zs <- seg.z[sort.seg]

  ##平均値をセグメントごとに発生させる
  (lower <- runif(col, 40, 80))   #変数の下限値
  (upper <- lower + runif(col, 30, 60))   #変数の上限値
  mu <- runif(col*seg[c], lower, upper)   #変数の平均値を発生
  (MU.s <- matrix(mu, nrow=seg[c], ncol=col, byrow=T))   #行列にまとめる
  
  covM.s <- array(0, dim=c(col, col, seg[c]))   #分散共分散行列の配列を初期化
  for(i in 1:seg[c]){  
    ##分散共分散行列をセグメントごとに発生させる
    #混合分布ごとの相関行列を作成
    
    cor <- corrM(col=col, lower=-0.5, upper=0.5)
    cov <- covmatrix(col=col, corM=cor, lower=min(MU.s[i, ])/1.5, upper=max(MU.s[i, ])/1.5)
    covM.s[, , i] <- cov$covariance
  
    ##変数を発生させる
    x.seg <- round(mvrnorm(length(seg.zs[seg.zs==i]), MU.s[i, ], covM.s[, , i]), 0)
    X.seg[[i]] <- x.seg 
  }
  
  Xlist[[c]] <-  X.seg 
  MU[[c]] <- MU.s
  covM[[c]] <- covM.s
  Z <- c(Z, seg.zs)
}

##リストからデータフレームに変換
X <- matrix(0, nrow=0, ncol=col)
for(c in 1:ca){
  X.cat <- do.call(rbind, Xlist[[c]])
  X <- rbind(X, X.cat)
}
Xz <- data.frame(seg=Z, cat=category, X=X)   #データを結合

##発生させた変数の要約
hist(Xz[Xz$cat==2, 3], breaks=25, col="grey", xlab="value", main="混合正規分布")   #分布をプロット
by(Xz[, 3:ncol(Xz)], Xz[, 1:2], function(x) round(colMeans(x), 2))   #カテゴリーセグメント別の平均
by(Xz[, 3:ncol(Xz)], Xz[, 1:2], function(x) round(var(x), 2))   #カテゴリーセグメント別の分散共分散行列
by(Xz[, 3:ncol(Xz)], Xz[, 1:2], function(x) round(cov2cor(var(x)), 2))   #カテゴリーセグメント別の相関行列
by(Xz[, 3:ncol(Xz)], Xz[, 1:2], function(x) summary(x))   #カテゴリーセグメント別の要約統計量
table(Xz[, 2], Xz[, 1])   #カテゴリーとセグメントのクロス集計

#カテゴリーとセグメントを無視した場合
round(colMeans(Xz[, 3:ncol(Xz)]), 2)
round(var(Xz[, 3:ncol(Xz)]), 2)
round(cor(Xz[, 3:ncol(Xz)]), 2)
summary(Xz[, 3:ncol(Xz)])


####EMアルゴリズムで混合正規分布モデル(Smart Shifter)を推定####
####EMアルゴリズムで用いる関数を定義####
##多変量正規分布の尤度関数
dmv <- function(x, mean.vec, S){
  LLo <- 1 / (sqrt((2 * pi)^nrow(S) * det(S))) *
         exp(-as.matrix(x - mean.vec) %*% solve(S) %*% t(x - mean.vec) / 2)
  return(LLo)
}

##観測データの対数尤度と潜在変数zの定義
LLobz <- function(X, seg, mean.M, S, r){
  LLind <- matrix(0, nrow(X), ncol=seg)   #対数尤度を格納する行列
  
  #混合多変量正規分布のセグメントごとの尤度を計算
  for(k in 1:seg){
    mean.vec <- mean.M[k, ]
    S_s <- S[[k]]
    Li <- apply(X, 1, function(x) dmv(x=t(x), mean.vec=mean.vec, S=S_s))   #多変量正規分布の尤度を計算
    LLi <- ifelse(Li==0, 10^-300, Li)
    LLind[, k] <- as.vector(LLi)
  }
  
  #対数尤度と潜在変数zの計算
  LLho <- matrix(r, nrow=nrow(X), ncol=seg, byrow=T) * LLind
  z <- LLho/matrix(apply(LLho, 1, sum), nrow=nrow(X), ncol=seg)   #zの計算 
  LLsum <- sum(log(apply(matrix(r, nrow=nrow(X), ncol=seg, byrow=T) * LLind, 1, sum)))   #観測データの対数尤度の和
  rval <- list(LLob=LLsum, z=z, LL=LLind, Li=Li)
  return(rval)
}

####カテゴリーごとにセグメント選択を含んだ混合多変量正規分布推定のEMアルゴリズム####
##カテゴリーでのベストセグメントのパラメータ推定値を格納する変数
aic.best <- c()
bic.best <- c()
seg.best <- c()
M.best <- list()
S.best <- list()
Z.best <- list()
r.best <- list()

for(bsc in 1:ca){
  #検討するセグメント数
  cat("カテゴリーは今", bsc, "ちかぁ\n", "認められないわぁ\n")
  seg_list <- c(2:6)  
  
  ##推定値を格納する変数を定義
  S.seg <- list()
  M.seg <- list()
  Z.seg <- list()
  r.seg <- list()
  AIC <- c()
  BIC <- c()
  
  ##セグメント数を変化させながら最適なセグメント決定
  for(bs in 1:length(seg_list)){
    s <- seg_list[bs] 
    print(s)
    
    ##初期値の設定
    ##kmeas法で初期値を設定
    XS <- Xz[Xz$cat==bsc, 3:ncol(Xz)]
    index.f <- kmeans(x=XS, s)$cluster
    
    #平均ベクトルと分散共分散行列の初期値をセグメントごとに逐次的に代入
    mean.M <- matrix(0, s, ncol(XS))
    S <- list()
    for(i in 1:s){
      mean.M[i, ] <- colMeans(XS[index.f==i, ])
      S[[i]] <- var(XS[index.f==i, ])
    }
   
    #混合率の初期値
    r <- as.numeric(table(index.f)/sum(table(index.f)))
    
    ##アルゴリズムの設定
    #対数尤度の初期化
    L <- LLobz(X=XS, seg=s, mean.M=mean.M, S=S, r=r)
    L1 <- L$LLob
    
    #更新ステータス
    dl <- 100   #EMステップでの対数尤度の差の初期化
    tol <- 0.1 
    iter <- 1
    max.iter <- 100
    Z.err <- c()
    
    ##EMアルゴリズムによる推定
    while(abs(dl) >= tol & iter <= max.iter){   #dlがtol以上の場合は繰り返す
      #Mステップの計算
      z <- L$z   #潜在変数zの出力
      
      #平均ベクトルと分散共分散行列を推定
      mean.M <- matrix(0, nrow=s, ncol=col)
      S <- list()
      for(js in 1:s){
        #平均ベクトルを推定
        mean.M[js, ] <- colSums(z[, js]*XS) / rep(nrow(XS)*r[js], col)
        
        #分散共分散行列を推定
        mean.v <- matrix(mean.M[js, ], nrow=nrow(XS), ncol=col, byrow=T)
        S[[js]] <- (t(z[, js]*as.matrix(XS) - z[, js]*mean.v) %*% 
                     (z[, js]*as.matrix(XS) - z[, js]*mean.v)) / sum(z[, js])
      }
      
      #混合率の推定
      r <- apply(L$z, 2, sum) / nrow(XS)
      
      ##Eステップの計算
      L <- try(LLobz(X=XS, seg=s, mean.M=mean.M, S=S, r=r), silent=TRUE)   #観測データの対数尤度を計算
      if(class(L) == "try-error") break   #エラー処理
      LL <- L$LLob   #観測データの対数尤度
      
      ##アルゴリズムの更新
      iter <- iter+1
      dl <- LL-LL1
      LL1 <- LL
      print(LL)
    }
    
    ##推定されたパラメータを格納
    #AICとBICの計算して格納
    AIC <- c(AIC, -2*LL + 2*(s*sum(1:col)+s*col))
    BIC <- c(BIC, -2*LL + log(nrow(XS))*(s*sum(1:col)+s*col))
    
    #パラメータを格納
    S.seg[[bs]] <- S
    M.seg[[bs]] <- mean.M
    if(class(try(L$z, silent=TRUE))=="try-error") {next} else {Z.seg[[bs]] <- L$z}   #errorの場合は次のセグメントへ
    r.seg[[bs]] <- r 
  }
  
  len <- length(Z.seg)   #エラーが起きていないセグメントを取得
  s.best <- which.min(AIC[1:len])   #AICで最適なセグメントを選択
  
  aic.best <- c(aic.best, AIC[s.best])
  bic.best <- c(bic.best, BIC[s.best])
  seg.best <- c(seg.best, (s.best+1))
  M.best[[bsc]] <- M.seg[[s.best]]
  S.best[[bsc]] <- S.seg[[s.best]]
  Z.best[[bsc]] <- Z.seg[[s.best]]
  r.best[[bsc]] <- r.seg[[s.best]]
}

####推定されたパラメータと真のセグメントの比較####
#カテゴリーの単純集計
table(Xz$cat) 
round(table(Xz$cat) / sum(table(Xz$cat)), 3)

##推定されたセグメントと真のセグメントを比較
#真のセグメント数を抽出
seg.t <- c()
for(i in 1:ca){
  seg.t <- c(seg.t, length(table(Xz[Xz$cat==i, 1])))
}
seg.best   #選択されたセグメント数
seg.t   #真のセグメント数

##平均ベクトルの比較
sapply(M.best, function(x) round(x, 2))   #推定された平均ベクトル
sapply(MU, function(x) round(x, 2))   #真の平均ベクトル

##分散共分散行列の比較
lapply(S.best[[1]], function(x) round(x, 2))   #推定された分散共分散行列
round(covM[[1]], 2)   #真の分散共分散行列

##潜在変数zの値
round(Z.best[[1]], 3)
round(Z.best[[2]], 3)
round(Z.best[[3]], 3)
round(Z.best[[4]], 3)
round(Z.best[[5]], 3)


