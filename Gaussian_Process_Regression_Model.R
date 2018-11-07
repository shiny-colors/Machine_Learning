######ガウス過程回帰モデル#####
options(warn=0)
library(MASS)
library(kernlab)
library(GPfit)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5698)

####データの発生####
#データの設定
d <- 1000   #学習データ
k <- 30   #入力変数数
w <- rpois(d, rgamma(d, 5.0, 0.7))

#テストデータのインデックス
index_learn <- 1:d1
index_test <- (d1+1):d

##データの生成
#時系列要素の変数を生成
#週成分を生成
week0 <- matrix(diag(7), nrow=d, ncol=7, byrow=T)
week <- week0[, -1]

#月成分を生成
m_days <- round(runif(ceiling(d/30), 30, 31))
month_vec <- rep(rep(1:12, round(length(m_days)/12))[1:length(m_days)], m_days)
month0 <- matrix(0, nrow=length(month_vec), ncol=12)
for(j in 1:12){
  index <- which(month_vec==j)
  month0[index, j] <- 1
}
month <- month0[1:d, -1]

#季節成分を生成
s_days <- rep(trunc(365/4), ceiling(d/(365/4))) + rbinom(ceiling(d/(365/4)), 1, 1/4)
season_vec <- rep(rep(1:4, length(s_days))[1:length(s_days)], s_days)
season0 <- matrix(0, nrow=length(season_vec), ncol=4)
for(j in 1:4){
  index <- which(season_vec==j)
  season0[index, j] <- 1
}
season <- season0[1:d, -1]


#多項分布から入力変数を生成
for(rp in 1:1000){
  Data0 <- matrix(0, nrow=d, ncol=k)
  for(j in 1:4){
    index <- which(season0[1:d, j]==1)
    alpha0 <- rep(0.3, k)
    theta <- extraDistr::rdirichlet(1, alpha0)   #多項分布のパラメータ
    Data0[index, ] <- rmnom(length(index), w[index], theta)
  }
  if(min(colSums(Data0)) >= 10) break
}

#データを結合
Data <- cbind(week, month, season, Data0)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")

#カーネル関数の生成
K <- Data %*% t(Data)

#ガウス過程からセグメントごとに応答変数を生成
sigma <- 2.0
y <- mvrnorm(1, rep(0, d), K + diag(sigma, d))   #ガウス過程から応答変数を生成
plot(1:d, y, type="l", xlab="time", ylab="y")


####EMアルゴリズムでガウス過程回帰モデルを推定####
##EMアルゴリズムの設定
iter <- 0
LL1 <- -10^100   #対数尤度の初期値
dl <- 100
tol <- 1

#初期値の設定
tau1 <- 1
tau2 <- 0.5

#データの設定
KK <- K %*% K


##EMアルゴリズムでパラメータを更新
while(abs(dl) >= tol){   #dlがtol以上なら繰り返す
  #Eステップで回帰係数を推定
  beta <- solve(KK + diag(tau1/tau2, d)) %*% t(K) %*% y
  delta <- diag(tau1, d) + tau2*KK
  
  #Mステップでハイパーパラメータを推定
  tau1_inv <- (sum(abs(beta^2)) + sum(diag(solve(delta)))) / d
  tau1 <- 1 / tau1_inv
  tau2_inv <- (sum((y - K %*% beta)^2) + sum(diag(KK %*% solve(delta)))) / d
  tau2 <- 1 / tau2_inv
  
  #周辺尤度の更新
  diag_tau2 <- diag(tau2_inv, d)
  tau1_KK <- tau1_inv * KK
  Lm <- -1/2*abs(diag_tau2 + tau1_KK) - 1/2*as.numeric((t(y) %*% solve(diag_tau2 + tau1_KK) %*% y))
  LL <- sum(Lm)
  
  ##EMアルゴリズムのパラメータの更新
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}


####推定結果の確認と要約####
##予測結果と実測値の比較
y_pred <- K %*% beta
plot(1:d, y, type="l", xlab="time", ylab="y", xlim=c(0, d), ylim=c(min(y), max(y)), main="実測値と予測値の時系列")
par(new=T)
plot(1:d, y_pred, type="l", col=2, xlab="", ylab="", xlim=c(0, d), ylim=c(min(y), max(y)))

##適合度を確認
sum((y - y_pred)^2)   #二乗誤差
sum(dnorm(y, y_pred, sd(y - y_pred), log=TRUE))   #対数尤度
