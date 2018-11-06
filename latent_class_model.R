#####潜在クラスモデル#####
####混合多項分布モデル####
####データの発生####
##モデルの設定
k <- 4   #セグメント数
n <- 500   #セグメントのサンプル数
N <- n * k   #総サンプル数
ns <- 30   #個人ごとの頻度
th <- 15   #セグメント別のパラメータ数

##確率ベクトルを定義
#セグメント1の出現率
x1 <- rnorm(th, 0, 1)
a1 <- rnorm(th, 0, 0.7)
(p1 <- round(exp(a1*x1) / sum(exp(a1*x1)), 3))

#セグメント2の出現率
x2 <- rnorm(th, 0, 1)
a2 <- rnorm(th, 0, 0.3)
(p2 <- round(exp(a2*x2) / sum(exp(a2*x2)), 3))

#セグメント3の出現率
x3 <- rnorm(th, 0, 1)
a3 <- rnorm(th, 0, 1.0)
(p3 <- round(exp(a3*x3) / sum(exp(a3*x3)), 3))

#セグメント4の出現率
x4 <- rnorm(th, 0, 1)
a4 <- rnorm(th, 0, 1.2)
(p4 <- round(exp(a4*x4) / sum(exp(a4*x4)), 3))

#真の確率
thetaT <- rbind(p1, p2, p3, p4)

##セグメントの確率に基づき頻度データを生成
#セグメントあたりの人数は500人、1人あたりの頻度は30個
#セグメント1のデータを発生
Y1 <- t(rmultinom(n, ns, p1))
dim(Y1)   #次元数

#セグメント2のデータを発生
Y2 <- t(rmultinom(n, ns, p2))
dim(Y2)   #次元数

#セグメント3のデータを発生
Y3 <- t(rmultinom(n, ns, p3))
dim(Y3)   #次元数

#セグメント4のデータを発生
Y4 <- t(rmultinom(n, ns, p4))
dim(Y4)   #次元数

#セグメントを表すベクトル
Seg <- rep(1:4, rep(n, 4))


##データを結合して結果を集計
Y <- as.data.frame(rbind(Y1, Y2, Y3, Y4))   #セグメントベクトルなしのデータ行列
Ys <- as.data.frame(cbind(Seg, Y))   #セグメントベクトルありのデータ行列
by(Ys[, 2:16], Ys[, 1], function(x) round(colMeans(x), 3))   #セグメント別の列ごとの発生頻度の平均
by(Ys[, 2:16], Ys[, 1], function(x) summary(x))   #セグメント別の集計
by(Ys[, 2:16], Ys[, 1], function(x) round(colSums(x)/sum(x), 3))   #セグメント別の発生率
dim(Y)   #次元数

####EMアルゴリズムで潜在クラスモデルを推定####
k <- 4   #セグメント数
n <- 500   #セグメントのサンプル数
N <- n * k   #総サンプル数
ns <- 30   #個人ごとの頻度
th <- 15   #パラメータ数

##観測データの対数尤度と潜在変数zを計算するための関数
LLobz <- function(theta, y, r, k){
  LLind <- matrix(0, nrow=nrow(y), ncol=k)
  for(i in 1:k){
    Li <- apply(y, 1, function(x) dmultinom(x, ns, theta[i, ]))   #多項分布の尤度を計算
    LLind[, i] <- Li
  }
  LLho <- matrix(r, nrow=nrow(y), ncol=k, byrow=T) * LLind   #観測データの尤度
  z <- LLho / matrix(apply(LLho, 1, sum), nrow=nrow(y), ncol=k)   #潜在変数zの計算
  LLosum <- sum(log(apply(matrix(r, nrow=nrow(y), ncol=k, byrow=T) * LLind, 1, sum)))   #観測データの対数尤度の計算
  rval <- list(LLob=LLosum, z=z, LL=LLind, Li=Li)
  return(rval)
}

#初期値の設定
iter <- 0
k <- 4   #セグメント数

##thetaの初期値の設定
#セグメント1の初期値
minmax <- colSums(Y)
hh1 <- runif(th, min(minmax), max(minmax))
theta1 <- hh1/sum(hh1)

#セグメント2の初期値
hh2 <- runif(th, min(minmax), max(minmax))
theta2 <- hh2/sum(hh2)

#セグメント3の初期値
hh3 <- runif(th, min(minmax), max(minmax))
theta3 <- hh3/sum(hh3)

#セグメント4の初期値
hh4 <- runif(th, min(minmax), max(minmax))
theta4 <- hh4/sum(hh4)

#thetaの初期値を結合
(theta <- rbind(theta1, theta2, theta3, theta4))

##混合率rの初期値
r <- c(0.4, 0.3, 0.2, 0.2)

#対数尤度の初期化
L <- LLobz(theta=theta, y=Y, r=r, k=k)
LL1 <- L$LLob
z <- L$z
round(z, 3)

#更新ステータス
dl <- 100   #EMステップでの対数尤度の差の初期値
tol <- 0.1  

##EMアルゴリズム
while(abs(dl) >= tol){   #dlがtol以上の場合は繰り返す
  #Eステップの計算
  z <- L$z   #潜在変数zの出力
  
  #Mステップの計算と最適化
  #thetaの推定
  theta <- matrix(0, nrow=k, ncol=th)
  for(j in 1:k){
    #完全データの対数尤度からthetaの推定量を計算
    thetaseg <- apply(matrix(z[, j], nrow=nrow(Y), ncol=th)*Y, 2, sum) / sum((z[, j])*matrix(ns, nrow=nrow(Y), ncol=1))
    theta[j, ] <- as.matrix(thetaseg)
  }
  #混合率を推定
  r <- apply(L$z, 2, sum) / nrow(Y)
  
  #観測データの対数尤度を計算
  L <- LLobz(theta=theta, y=Y, r=r, k=k)
  LL <- L$LLob   #観測データの対数尤度
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####潜在クラスモデルの推定結果####
round(theta, 3)   #thetaの推定量
round(thetaT, 3)   #thetaの真の値
round(r, 3)   #混合率の推定量
round(z, 3)   #個人別のセグメントへの所属確率

L$LLob   #観測データの対数尤度
-2*(L$LLob) + k*nrow(theta)   #AIC

