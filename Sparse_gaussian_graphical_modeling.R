#####スパースガウシアングラフィカルモデリング#####
library(MASS)
library(lars)
library(glmnet)
library(glasso)
library(matrixStats)
library(Matrix)
library(bayesm)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(igraph)

####多変量正規分布の乱数を発生させる関数を定義####
#任意の相関行列を作る関数を定義
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #新しい相関行列の定義と対角成分を1にする
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
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
##データの設定
N <- 10000   #サンプル数
k <- 20   #変数数
mu <- rep(0, k)   #平均ベクトル

##多変量正規分布からデータを生成
Cor <- corrM(k, -0.6, 1.0, 0.1, 1.0)
Data <- mvrnorm(N, mu, Cor)
S <- cor(Data)   #観測相関行列


#####グラフィカルlassoでガウシアングラフィカルモデルを推定####
##アルゴリズムの設定
tol <- 1
lambda <- 0.04   #正則化パラメータ
diff <- 100

#初期値の設定
Sigmat <- Sigma <- S + diag(lambda, k)
Ohm <- solve(Sigma)
dl <- sum(abs(Sigma[upper.tri(Sigma)]))

#インデックスを作成
index1 <- matrix(0, nrow=k, ncol=k)
index2 <- matrix(0, nrow=k, ncol=k)
for(i in 1:k){
  index1[i, ] <- c((1:k)[-i], i)
  index2[i, i] <- k
  index2[i, -i] <- 1:(k-1)
}

##グラフィカルlassoでパラメータを更新
while(diff > tol){
  for(i in 1:k){
  
    #入力変数の入れ替え
    Ohm_tilde <- Ohm[index1[i, ], index1[i, ]]^1/2
    Sigma_tilde <- Sigma[index1[i, ], index1[i, ]]
    S_tilde <- S[index1[i, ], index1[i, ]]
    
    #lasso回帰でパラメータベクトルを推定
    y <- as.numeric((solve(Sigma_tilde[-k, -k])^1/2) %*% S_tilde[-k, k])
    X <- Sigma_tilde[-k, -k]^1/2
    res <- glmnet(X, y, family="gaussian", lambda=lambda, intercept=FALSE, standardize=TRUE)
    
    #分散共分散行列を推定
    beta <- as.numeric(Sigma[-i, -i] %*%  res$beta)
    Sigma[-i, i] <- Sigma[i, -i] <- beta
    omega <- as.numeric(1 / (Sigma_tilde[k, k] - t(beta) %*% res$beta))
    Ohm[-i, i] <- Ohm[i, -i] <- -omega * beta
  }
  
  #収束判定
  dl1 <- sum(abs(Sigma[upper.tri(Sigma)]))
  diff <- abs(dl1 - dl)
  dl <- dl1
  print(diff)
}

####相関構造を可視化####
diag(Sigma) <- 0
Sigma <- ifelse(Sigma > 0, 1, ifelse(Sigma < 0, -1, 0))
g <- graph.adjacency(Sigma, mode="undirected")
plot(g, vertex.size=20, vertex.shape="rectangle", vertex.color="#FFFF99")


##関数で推定
res <- glasso(S, lambda)
round(solve(res$wi), 2)

