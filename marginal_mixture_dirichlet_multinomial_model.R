#####周辺化多項ディクレリ混合モデル#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####データの発生####
##モデルの設定
k <- 5   #セグメント数
hh <- 500   #セグメントのサンプル数
n <- hh * k   #総サンプル数
w <- 200   #セグメント別のパラメータ数

#セグメントの設定
seg.z <- rep(1:k, rep(hh, k))

##個人ごとの頻度を設定
freq <- round(exp(rnorm(n, 4.75, 0.18)), 0)
summary(freq)
hist(freq, col="grey", main="ユーザーごとの閲覧数")

##セグメントごとに確率ベクトルを定義
#ディクレリ分布より確率係数を設定
x <- rdirichlet(k, rep(0.1, w))
round(x, 3)   #発生させた確率を確認

##発生させた確率にもとづき頻度データを生成
#多項分布より閲覧履歴を生成
Y <- matrix(0, nrow=n, ncol=w)
for(i in 1:n){
  Y[i, ] <- t(rmultinom(1, freq[i], x[seg.z[i], ]))
}
matrix(as.numeric(unlist(by(Y, seg.z, colSums))), nrow=k, ncol=w, byrow=T)   #セグメントごとの出現数を計算


####周辺化ギブスサンプリングで多項ディクレリ分布のセグメントを生成####
##MCMCアルゴリズムの設定
R <- 5000
keep <- 2
sbeta <- 1.5

##ハイパーパラメータの設定
alpha <- 5
beta <- 5
beta.array <- array(beta, dim=c(k, w, n))


####周辺化ギブスサンプリングでセグメントを生成####
##zの初期値を設定
#初期値依存するので、良い初期値が得られるまで反復させる
for(t in 1:1000){
  print(t)
  #ユーザーごとのセグメント割当を行列形式に変更
  Z <- t(rmultinom(n, 1, rep(1/k, k)))
  z <- Z %*% 1:k
  
  zu <- matrix(0, nrow=n, ncol=k)
  za <- matrix(0, nrow=n, ncol=k)
  zw <- array(0, dim=c(k, w, n))
  
  for(i in 1:n){
    zu[i, z[i]] <- 1 
    za[i, z[i]] <- freq[i]
    zw[z[i], , i] <- Y[i, ]
  }
  
  ##データとアルゴリズム仮定の格納用パラメータを設定
  Nd <- matrix(freq, nrow=n, ncol=k)   #ユーザーごとの頻度
  
  #ユーザーおよびアイテムごとの頻度
  Ndv <- array(0, dim=c(k, w, n))   
  for(i in 1:n){
    Ndv[, , i] <- matrix(Y[i, ], nrow=k, ncol=w, byrow=T)
  }
  Fdv <- ifelse(Ndv > 0, 1, 0)
  
  #ユーザーごとにアイテムの出現がある列を取得
  index.n <- list()
  for(i in 1:n){
    index.n[[i]] <- subset(1:length(Y[i, ]), Y[i, ] > 0)
  }
  
  #ユーザーごとの頻度行列
  freqM <- matrix(freq, nrow=n, ncol=k)
  
  #頻度の多次元配列
  Y_array <- array(0, dim=c(k, w, n))
  for(i in 1:n){
    Y_array[, , i] <- matrix(Y[i, ], nrow=k, ncol=w, byrow=T)
  }
  
  #ポリア分布のパラメメータ用
  lgamma3 <- matrix(0, nrow=n, ncol=k)
  lgamma4 <- matrix(0, nrow=n, ncol=k)
  
  ##サンプリング結果の格納用
  Seg_Z <- matrix(0, nrow=R/keep, ncol=n)
  Prob <- array(0, dim=c(n, k, R/keep))
  
  ##統計量の初期値を計算
  #ユーザーのセグメント割当を解除
  Zu <- matrix(as.numeric(table(z)), nrow=n, ncol=k, byrow=T)
  Zul <- Zu - zu
  
  #総頻度のセグメント割当を解除
  Za <- matrix(as.numeric(by(Y, z, sum)), nrow=n, ncol=k, byrow=T)
  Zal <- Za - za
  
  #wのセグメント割当を解除
  Zw <- array(matrix(unlist(by(Y, z, colSums)), nrow=k, ncol=w, byrow=T), dim=c(k, w, n))
  Zwl <- Zw - zw
  
  
  ##ここからmcmcサンプリングを実行
  for(rp in 1:R){
  
    ##ギブスサンプリングの確率の更新式を計算
    #第1因子を計算
    D_alpha <- log(Zul + alpha)
    
    #第2因子を計算
    lgamma1 <- lgamma(Zal + beta*w)
    lgamma2 <- lgamma(Zal + Nd + beta*w)
    
    #第3因子を計算
    lg4 <- Zwl + beta.array*Fdv
    lg3 <- lg4 + Ndv 
    
    for(i in 1:n){
      lgamma3[i, ] <- rowSums(lgamma(lg3[, index.n[[i]], i]))
      lgamma4[i, ] <- rowSums(lgamma(lg4[, index.n[[i]], i]))
    }
    
    ##確率の計算とセグメントの生成
    #セグメント割当確率を計算
    #尤度が桁落ちしないように定数を加える
    lgamma_c <- D_alpha + lgamma1 - lgamma2 + lgamma3 - lgamma4
    lgamma_mean <- matrix(apply(lgamma_c, 1, mean), nrow=n, ncol=k)
    lgamma <- exp(lgamma_c - lgamma_mean)
    lgamma[is.infinite(lgamma) & lgamma > 0] <- 10^300
    lgamma[is.infinite(lgamma) & lgamma < 0] <- -10^300
    
    #確率の計算
    Pr <- lgamma / rowSums(lgamma)
  
    #多項分布よりセグメントを生成
    Z <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
    z <- Z %*% 1:k
    
    if(length(table(z)) < k) {break}   #セグメントが縮退したらbreakする
    
    ##統計量を更新
    #セグメント割当を解除
    #ユーザーのセグメント割当を解除
    Zu <- matrix(as.numeric(table(z)), nrow=n, ncol=k, byrow=T)
    zu <- Z
    Zul <- Zu - zu
    
    #総頻度のセグメント割当を解除
    Za <- matrix(as.numeric(by(Y, z, sum)), nrow=n, ncol=k, byrow=T)
    za <- Z * freqM   
    Zal <- Za - za
    
    #wのセグメント割当を解除
    Zw <- array(matrix(unlist(by(Y, z, colSums)), nrow=k, ncol=w, byrow=T), dim=c(k, w, n))
    Z_array <- array(as.numeric(t(Z[rep(1:n, rep(w, n)), ])), dim=c(k, w, n))   
    zw <- Z_array * Y_array
    Zwl <- Zw - zw
    
    ##パラメータを格納
    if(rp%%keep==0){
      print(rp)
      mkeep <- rp/keep
      Seg_Z[mkeep, ] <- z
      Prob[, , mkeep] <- Pr
      print(table(z))
    }
  }
  if(rp==R) {break}   #設定したサンプリング数に達していれば、mcmcループを終了する
}

####推定結果の要約と適合度の確認####
burnin <- 2000/keep   #バーンイン期間

#セグメントのサンプリング結果
Seg_Z[burnin:(R/keep), c(1:4, 501:504, 1001:1004, 1501:1504, 2001:2004)]

##確率のサンプリング結果
round(colMeans(t(Prob[1, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[501, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[1001, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[1501, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[2001, , burnin:(R/keep)])), 3)
