#####Tree Structured Mixture Multinomial Model#####
options(warn=0)
library(stringr)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
#set.seed(2506787)

####データの発生####
##すべての単語が出現するまでデータの生成を続ける
for(rp in 1:1000){
  print(rp) 
  
  ##データの設定
  #木構造のトピックの設定
  for(iter in 1:1000){
    m <- 3   #木の最大深さ
    k1 <- 3   #1階層目のトピック数
    k2 <- rtpois(k1, a=1, b=5, 2.5)   #2階層目のトピック数
    k3 <- c()
    for(j in 1:k1){
      k3 <- c(k3, rtpois(k2[j], a=0, b=4, 2.25))   #3階層目のトピック数
    }
    index_k3 <- cbind(c(1, cumsum(k2)[-length(k2)]+1), cumsum(k2))
    k <- 1 + sum(c(k1, k2, k3))   #総トピック数
    max_k <- max(c(k1, k2, k3))
    if(k >= 30 & k < 40) break 
  }
  
  #文書の設定
  d <- 3500   #文書数
  w <- rpois(d, rgamma(d, 75, 0.5))   #単語数
  f <- sum(w)   #総単語数
  v0 <- 250
  v1 <- rep(250, k1)
  v <- sum(c(v0, v1))   #総語彙数
  index_v0 <- 1:v0
  index_v1 <- cbind(c(v0+1, (v0+cumsum(v1)[-k1])+1), (v0+cumsum(v1)))
  
  #IDを設定
  d_id <- rep(1:d, w)
  t_id <- c()
  for(i in 1:d){
    t_id <- c(t_id, 1:w[i])
  }
  
  ##パラメータを生成
  #単語分布の事前分布
  alpha1 <- rep(0.01, v)
  alpha21 <- alpha22 <- alpha23 <- rep(0.0005, v)
  alpha1[index_v0] <- 5.0
  alpha21[index_v1[1, 1]:index_v1[1, 2]] <- 0.1
  alpha22[index_v1[2, 1]:index_v1[2, 2]] <- 0.1
  alpha23[index_v1[3, 1]:index_v1[3, 2]] <- 0.1
  alpha2 <- rbind(alpha21, alpha22, alpha23)
  
  #トピック割当の事前分布
  beta1 <- c(30.0, 20.0)   #木の深さの事前分布
  beta2 <- rep(5.0, max(c(k1, k2, k3)))   #分岐の事前分布

  ##木構造のパラメータを生成
  #ベータ分布から停止確率を生成
  gamma1 <- gammat1 <- 0.8
  gamma2 <- gammat2 <- rbeta(k1, beta1[1], beta1[2])
  gamma3 <- list()
  for(j in 1:k1){
    gamma3[[j]] <- rbeta(k2[j], beta1[1], beta1[2])
  }
  gammat3 <- gamma3
  
  #ディリクリ分布から木構造のノード選択確率を生成
  theta1 <- thetat1 <- as.numeric(extraDistr::rdirichlet(1, rep(20.0, k1)))
  theta2 <- theta3 <- list()
  for(i in 1:k1){
    theta_list <- list()
    theta2[[i]] <- as.numeric(extraDistr::rdirichlet(1, beta2[1:k2[i]]))
    
    for(j in 1:k2[i]){
      x <- k3[(index_k3[i, 1]:index_k3[i, 2])][j]
      if(x==1){
        theta_list[[j]] <- 1
      } else {
        theta_list[[j]] <- as.numeric(extraDistr::rdirichlet(1, beta2[1:x]))
      }
    }
    theta3[[i]] <- theta_list
  }
  thetat2 <- theta2; thetat3 <- theta3
  
  ##単語分布を生成
  phi0 <- phit0 <- as.numeric(extraDistr::rdirichlet(1, alpha1))
  phi1 <- phit1 <- extraDistr::rdirichlet(k1, alpha2)
  phi2 <- phi3 <- list()
  
  for(i in 1:k1){
    phi_list <- list()
    phi2[[i]] <- extraDistr::rdirichlet(k2[i], alpha2[i, ])
    
    for(j in 1:k2[i]){
      x <- k3[(index_k3[i, 1]:index_k3[i, 2])][j]
      phi_list[[j]] <- extraDistr::rdirichlet(x, alpha2[i, ])
    }
    phi3[[i]] <- phi_list
  }
  phit2 <- phi2; phit3 <- phi3
  
  ##モデルに基づきデータを生成
  #データの格納用配列
  Z1_list <- Z2_list <- word_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    ##木構造のノードと通過を生成
    #初期ノードの通過を生成
    z11_vec <- rbinom(w[i], 1, gamma1)   
    
    #1階層目のトピックを生成
    z21 <- rmnom(w[i], 1, theta1) * z11_vec
    z21_vec <- as.numeric(z21 %*% 1:k1)
    
    #1階層目のノードの通過を生成
    z12_vec <- rep(0, w[i])
    z12_vec[z11_vec==1] <- rbinom(sum(z11_vec), 1, gamma2[z21_vec])
    
    #2階層目のトピックを生成
    n <- sum(z12_vec)
    x21 <- z21_vec[z12_vec==1]
    theta_par <- matrix(0, nrow=n, ncol=max_k)
    for(j in 1:n){
      theta_par[j, 1:k2[x21[j]]] <- theta2[[x21[j]]]
    }
    z22 <- matrix(0, nrow=w[i], ncol=max_k)
    z22[z12_vec==1, ] <- rmnom(n, 1, theta_par)   #多項分布からトピックを生成
    z22_vec <- as.numeric(z22 %*% 1:max_k)
    
    #2階層目のノードの通過を生成
    gamma_par <- rep(0, n)
    z22_vec1 <- z22_vec[z22_vec > 0]
    for(j in 1:n){
      gamma_par[j] <- gamma3[[x21[j]]][z22_vec1[j]]
    }
    z13_vec <- rep(0, w[i])
    z13_vec[z12_vec==1] <- rbinom(n, 1, gamma_par)
    
    #3階層目のトピックを生成
    n <- sum(z13_vec)
    x22 <- cbind(x21[z13_vec[z12_vec==1]==1], z22_vec[z13_vec==1])
    theta_par <- matrix(0, nrow=n, ncol=max_k)
    
    for(j in 1:n){
      k0 <- k3[(index_k3[x22[j, 1], 1]:index_k3[x22[j, 1], 2])[x22[j, 2]]]
      if(k0==1){
        theta_par[j, 1] <- 1
      } else {
        theta_par[j, 1:k0] <- theta3[[x22[j, 1]]][[x22[j, 2]]]
      }
    }
    z23 <- matrix(0, nrow=w[i], ncol=max_k)
    z23[z13_vec==1, ] <- rmnom(n, 1, theta_par)   #多項分布からトピックを生成
    z23_vec <- as.numeric(z23 %*% 1:max_k)
    
    #トピックを結合
    z2 <- cbind(z0=1, z21_vec, z22_vec, z23_vec)
    depth <- rowSums(z2 > 0)
    
    
    ##単語を生成
    #木構造に応じて単語分布のパラメータを格納
    phi <- matrix(0, nrow=w[i], ncol=v)
    phi[depth==1, ] <- matrix(phi0, nrow=sum(1-z11_vec), ncol=v, byrow=T)
    phi[depth==2, ] <- phi1[z2[depth==2, 2], ]
    
    for(j in 1:w[i]){
      if(depth[j] < 3) next
      if(depth[j]==3){
        phi[j, ]  <- phi2[[z2[j, 2]]][z2[j, 3], ]
      }
      if(depth[j]==4){
        phi[j, ] <- phi3[[z2[j, 2]]][[z2[j, 3]]][z2[j, 4], ]
      }
    }
    
    #多項分布から単語を生成
    word <- rmnom(w[i], 1, phi)
    word_vec <- as.numeric(word %*% 1:v)
    
    ##データの格納
    Z1_list[[i]] <- cbind(z11_vec, z12_vec, z13_vec)
    Z2_list[[i]] <- z2
    word_list[[i]] <- word_vec
    WX[i, ] <- colSums(word)
  }
  if(min(colSums(WX)) > 0) break
}

#リストを変換
wd <- unlist(word_list)
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
sparse_data <- sparseMatrix(i=1:f, wd, x=rep(1, f), dims=c(f, v))
sparse_data_T <- t(sparse_data)


####マルコフ連鎖モンテカルロ法でTree Structured Mixture Multinomial Modelを推定####
##アルゴリズムの設定
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##事前分布の設定
alpha011 <- 1000.0
alpha012 <- 20.0
alpha02 <- 0.01
alpha03 <- 0.01
s0 <- 0.5
v0 <- 0.5

##パラメータの真値
gamma1 <- gammat1
gamma2 <- gammat2
gamma3 <- gammat3
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3
phi0 <- phit0
phi1 <- phit1
phi2 <- phit2
phi3 <- phit3

##パラメータの初期値
#ベータ分布から停止確率を生成
gamma1 <- 0.7
gamma2 <- rbeta(k1, 20, 15)
gamma3 <- list()
for(j in 1:k1){
  gamma3[[j]] <- rbeta(k2[j],20, 15)
}

#ディリクリ分布から木構造のノード選択確率を生成
theta1 <- as.numeric(extraDistr::rdirichlet(1, rep(1/k1, k1)))
theta2 <- theta3 <- list()
for(i in 1:k1){
  theta_list <- list()
  theta2[[i]] <- as.numeric(extraDistr::rdirichlet(1, rep(1/k2[i], k2[i])))

  for(j in 1:k2[i]){
    x <- k3[(index_k3[i, 1]:index_k3[i, 2])][j]
    if(x==1){
      theta_list[[j]] <- 1
    } else {
      theta_list[[j]] <- as.numeric(extraDistr::rdirichlet(1, rep(1/length(1:x), length(1:x))))
    }
  }
  theta3[[i]] <- theta_list
}

##単語分布を生成
phi0 <- as.numeric(extraDistr::rdirichlet(1, rep(10.0, v)))
phi1 <- extraDistr::rdirichlet(k1, rep(10.0, v))
phi2 <- phi3 <- list()

for(i in 1:k1){
  phi_list <- list()
  phi2[[i]] <- extraDistr::rdirichlet(k2[i], rep(10.0, v))
  
  for(j in 1:k2[i]){
    x <- k3[(index_k3[i, 1]:index_k3[i, 2])][j]
    phi_list[[j]] <- extraDistr::rdirichlet(x, rep(10.0, v))
  }
  phi3[[i]] <- phi_list
}

##パラメータの格納用配列
PHI0 <- matrix(0, nrow=R/keep, ncol=v)
PHI1 <- array(0, dim=c(k1, v, R/keep))
PHI2 <- array(0, dim=c(sum(k2), v, R/keep))
PHI3 <- array(0, dim=c(sum(k3), v, R/keep))
GAMMA1 <- rep(0, R/keep)
GAMMA2 <- matrix(0, nrow=R/keep, ncol=k1)
GAMMA3 <- matrix(0, nrow=R/keep, ncol=sum(k2))
THETA1 <- matrix(0, nrow=R/keep, ncol=k1)
THETA2 <- matrix(0, nrow=R/keep, ncol=sum(k2))
THETA3 <- matrix(0, nrow=R/keep, ncol=sum(k3))
SEG21 <- matrix(0, nrow=f, ncol=k1)
SEG22 <- matrix(0, nrow=f, ncol=sum(k2))
SEG23 <- matrix(0, nrow=f, ncol=sum(k3))

##インデックスを作成
doc_list <- doc_vec <- list()
wd_list <- wd_vec <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}
for(j in 1:v){
  wd_list[[j]] <- which(wd==j)
  wd_vec[[j]] <- rep(1, length(wd_list[[j]]))
}

##対数尤度の基準値
LLst <- sum(sparse_data %*% log(colSums(sparse_data) / f))


####ギブスサンプリングでパラメータをサンプリング####
for(rp in 1:R){
  
  ##木構造のノードごとの期待尤度の設定
  Lho0 <- phi0[wd]
  Lho1 <- matrix(theta1, nrow=f, ncol=k1, byrow=T) * t(phi1)[wd, ]
  Lis1 <- as.numeric((matrix(1-gamma2, nrow=f, ncol=k1, byrow=T) * Lho1) %*% rep(1, k1))
  Lho2 <- list(); Li2 <- matrix(0, nrow=f, ncol=k1)
  Lho3 <- Li3 <- list()
  
  for(i in 1:k1){
    #2階層目の期待尤度
    Lho2[[i]] <- matrix(theta2[[i]], nrow=f, ncol=k2[i], byrow=T) * t(phi2[[i]])[wd, ]
    r <- matrix((1-gamma3[[i]]), nrow=f, ncol=length(gamma3[[i]]), byrow=T)
    Li2[, i] <- as.numeric((r * Lho2[[i]]) %*% rep(1, k2[i]))
    
    #3階層目の期待尤度
    index <- (index_k3[i, 1]:index_k3[i, 2])
    a <- k3[index]
    Lho3_list <- list(); Li3_data <- matrix(0, nrow=f, ncol=k2[i])
    for(j in 1:k2[i]){
      Lho3_list[[j]] <- matrix(theta3[[i]][[j]], nrow=f, ncol=a[j], byrow=T) * t(phi3[[i]][[j]])[wd, ]
      Li3_data[, j] <- as.numeric(Lho3_list[[j]] %*% rep(1, a[j]))
    }
    Lho3[[i]] <- Lho3_list
    Li3[[i]] <- Li3_data
  }
  
  ##通過変数を生成
  #通過確率の尤度を設定
  Lis21 <- as.numeric((matrix(gamma2*theta1, nrow=f, ncol=k1, byrow=T) * Li2) %*% rep(1, k1))
  Lis31 <- matrix(0, nrow=f, ncol=k1)
  for(j in 1:k1){
    Lis31[, j] <- as.numeric((matrix(gamma3[[j]]*theta2[[j]], nrow=f, k2[j], byrow=T) * Li3[[j]]) %*% rep(1, k2[j]))
  }
  Lis32 <- as.numeric((matrix(gamma2*theta1, nrow=f, ncol=k1, byrow=T) * Lis31) %*% rep(1, k1))
  
  
  #1階層目の通過変数を生成
  Lis_par1 <- gamma1 * (Lis1 + Lis21 + Lis32)
  Lis_par0 <- (1-gamma1) * Lho0
  s_prob1 <- Lis_par1 / (Lis_par1 + Lis_par0)   #通過確率
  Zi11 <- rbinom(f, 1, s_prob1)   #ベルヌーイ分布から通過変数を生成
  index_z11 <- which(Zi11==1)
  
  #2階層目の通過変数を生成
  Zi12 <- s_prob2 <- rep(0, f)
  Lis_par1 <- (Lis21 + Lis32)[index_z11]
  Lis_par0 <- Lis1[index_z11]
  s_prob2[index_z11] <- Lis_par1 / (Lis_par1 + Lis_par0)   #通過確率
  Zi12[index_z11] <- rbinom(length(index_z11), 1, s_prob2[index_z11])   #ベルヌーイ分布から通過変数を生成
  index_z12 <- which(Zi12==1)
  
  #3階層目の通過変数を生成
  Zi13 <- s_prob3 <- rep(0, f)
  Lis_par1 <- (as.numeric(Lis31 %*% rep(1, k1)))[index_z12]
  Lis_par0 <- (as.numeric(Li2 %*% rep(1, k1)))[index_z12]
  s_prob3[index_z12] <- Lis_par1 / (Lis_par1 + Lis_par0)   #通過確率
  Zi13[index_z12] <- rbinom(length(index_z12), 1, s_prob3[index_z12])   #ベルヌーイ分布から通過変数を生成
  index_z13 <- which(Zi13==1)
  
  
  ##通過変数を条件づけたトピックをサンプリング
  ##1階層目のトピックを生成
  #ノード割当を設定
  topic_prob1 <- matrix(0, nrow=f, ncol=k1)
  index_d10 <- index_z11[Zi12[index_z11]==0]
  index_d11 <- index_z11[Zi12[index_z11]==1]
  
  #トピックの潜在割当確率
  r <- matrix(gamma2*theta1, nrow=f, ncol=k1, byrow=T)[index_d11, ]
  par <- (r*Li2[index_d11, ]) + (r*Lis31[index_d11, ])
  topic_prob1[index_d10, ] <- Lho1[index_d10, ] / as.numeric((Lho1[index_d10, ]) %*% rep(1, k1))   #トピックの割当確率
  topic_prob1[index_d11, ] <- par / rowSums(par)
  
  #多項分布からトピックを生成
  Zi21 <- matrix(0, nrow=f, ncol=k1)
  Zi21[index_d10, ] <- rmnom(length(index_d10), 1, topic_prob1[index_d10, ])   
  Zi21[index_d11, ] <- rmnom(length(index_d11), 1, topic_prob1[index_d11, ])
  
  ##2階層目のトピックを生成
  Zi22_list <- list(); Zi23_list <- list()
  for(i in 1:k1){
    #ノード割当を設定
    topic_prob2 <- matrix(0, nrow=f, ncol=k2[i])
    index_d20 <- index_z12[Zi13[index_z12]==0 & Zi21[index_z12, i]==1]
    index_d21 <- index_z12[Zi13[index_z12]==1 & Zi21[index_z12, i]==1]
    
    #トピックの潜在割当確率
    par20 <- (Lho2[[i]])[index_d20, ]
    par21 <- (matrix(gamma3[[i]]*theta2[[i]], nrow=f, k2[i], byrow=T) * Li3[[i]])[index_d21, ]
    topic_prob2[index_d20, ] <- par20 / rowSums(par20)   #トピックの割当確率
    topic_prob2[index_d21, ] <- par21 / rowSums(par21)
    
    #多項分布からトピックを生成
    Zi22 <- matrix(0, nrow=f, ncol=k2[i])
    Zi22[index_d20, ] <- rmnom(length(index_d20), 1, topic_prob2[index_d20, ])
    Zi22[index_d20, ][is.na(Zi22[index_d20, ])] <- 0
    Zi22[index_d21, ] <- rmnom(length(index_d21), 1, topic_prob2[index_d21, ])
    Zi22[index_d21, ][is.na(Zi22[index_d21, ])] <- 0
    
    ##3層目のトピックを生成
    Zi23_list0 <- list()
    for(j in 1:k2[i]){
      
      #ノード割当を設定
      index_d3 <-index_z13[Zi22[index_z13, j]==1]
      if(ncol(Lho3[[i]][[j]])==1){
        Zi23 <- matrix(0, nrow=f, ncol=1)
        Zi23[index_d3, ] <- 1
        Zi23_list0[[j]] <- Zi23
        next
      }
      #トピックの潜在割当確率
      par3 <- Lho3[[i]][[j]][index_d3, , drop=FALSE]
      topic_prob3 <- matrix(0, nrow=f, ncol=ncol(par3))
      topic_prob3[index_d3, ] <- par3 / rowSums(par3)   #トピックの割当確率
      
      #多項分布からトピックを生成
      Zi23 <- matrix(0, nrow=f, ncol=ncol(par3))
      Zi23[index_d3, ] <- rmnom(length(index_d3), 1, topic_prob3[index_d3, ])
      Zi23_list0[[j]] <- Zi23
    }
    Zi22_list[[i]] <- Zi22
    Zi23_list[[i]] <- Zi23_list0
  }
  
  ##パラメータをサンプリング
  ##0階層目のパラメータをサンプリング
  #通過率をサンプリング
  n1 <- sum(Zi11)
  s1 <- n1 + s0; v1 <- f - n1 + v0   #ベータ分布のパラメータ
  gamma1 <- rbeta(1, s1, v1)   #パラメータをサンプリング
  
  #単語分布をサンプリング
  vsum0 <- as.numeric(sparse_data_T %*% (1-Zi11)) + alpha02   #ディリクレ分布のパラメータ
  phi0 <- as.numeric(extraDistr::rdirichlet(1, vsum0))   #パラメータをサンプリング
  
  
  ##1階層目のパラメータをサンプリング
  #通過率をサンプリング
  n21 <- colSums(Zi21); n22 <- as.numeric(Zi12 %*% Zi21)
  s2 <- n22 + s0; v2 <- n21 - n22 + v0   #ベータ分布のパラメータ
  gamma2 <- rbeta(k1, s2, v2)   #パラメータをサンプリング
  
  #混合率をサンプリング
  wsum1 <- colSums(Zi21) + alpha011   #ディリクレ分布のパラメータ 
  theta1 <- as.numeric(extraDistr::rdirichlet(1, wsum1))   #パラメータをサンプリング
  
  #単語分布をサンプリング
  vsum1 <- t(sparse_data_T %*% (Zi21 * (1-Zi12)*Zi11)) + alpha02   #ディリクレ分布のパラメータ
  phi1 <- extraDistr::rdirichlet(k1, vsum1)
  
  
  ##2階層目のパラメータをサンプリング
  Zi_path <- Zi11 * Zi12   #パス変数
  
  for(i in 1:k1){
    #通過率をサンプリング
    n21 <- colSums(Zi22_list[[i]]); n22 <- as.numeric(Zi13 %*% Zi22_list[[i]])
    s3 <- n22 + s0; v3 <- n21 - n22 + v0   #ベータ分布のパラメータ
    gamma3[[i]] <- rbeta(length(n21), s3, v3)
    gammat3[[i]]
    
    #混合率をサンプリング
    wsum2 <- colSums(Zi22_list[[i]]) + alpha012   #ディリクレ分布のパラメータ
    theta2[[i]] <- as.numeric(extraDistr::rdirichlet(1, wsum2))   #パラメータをサンプリング
    
    #単語分布をサンプリング
    vsum2 <- t(sparse_data_T %*% (Zi22_list[[i]] * (1-Zi13)*Zi_path)) + alpha02
    phi2[[i]] <- extraDistr::rdirichlet(k2[i], vsum2)
    
    ##3階層目のパラメータをサンプリング
    for(j in 1:k2[i]){
      #混合率をサンプリング
      if(ncol(Zi23_list[[i]][[j]])==1){
        theta3[[i]][[j]] <- 1
      } else {
      wsum3 <- colSums(Zi23_list[[i]][[j]]) + alpha012   #ディリクレ分布のパラメータ
      theta3[[i]][[j]] <- extraDistr::rdirichlet(1, wsum3)   #パラメータをサンプリング
      }
      
      #単語分布をサンプリング
      vsum3 <- t(sparse_data_T %*% (Zi23_list[[i]][[j]] * Zi13*Zi_path)) + alpha03
      phi3[[i]][[j]] <- extraDistr::rdirichlet(nrow(vsum3), vsum3)   #パラメータをサンプリング
      }
  }
  
  ##パラメータの格納とサンプリング結果の表示
  #サンプリングされたパラメータを格納
  if(rp%%keep==0){
    #サンプリング結果の格納
    mkeep <- rp/keep
    PHI0[mkeep, ] <- phi0
    PHI1[, , mkeep] <- phi1
    PHI2[, , mkeep] <- do.call(rbind, phi2)
    PHI3[, , mkeep] <- rbind(do.call(rbind, phi3[[1]]), do.call(rbind, phi3[[2]]), do.call(rbind, phi3[[3]]))
    GAMMA1[mkeep] <- gamma1
    GAMMA2[mkeep, ] <- gamma2
    GAMMA3[mkeep, ] <- unlist(gamma3)
    THETA1[mkeep, ] <- theta1
    THETA2[mkeep, ] <- unlist(theta2)
    THETA3[mkeep, ] <- unlist(theta3)
  }
  
  #トピック割当はバーンイン期間を超えたら格納する
  if(rp%%keep==0 & rp >= burnin){
    SEG21 <- SEG21 + Zi21
    SEG22 <- SEG22 + do.call(cbind, Zi22_list)
    SEG23 <- SEG23 + cbind(do.call(cbind, Zi23_list[[1]]), do.call(cbind, Zi23_list[[2]]), do.call(cbind, Zi23_list[[3]]))
  }

  ##対数尤度の計算とサンプリング結果の確認
  if(rp%%disp==0){
    ##対数尤度を計算
    #データの設定
    Zi1 <- cbind(Zi11, Zi12, Zi13)
    index <- which(Zi1[, 1]==1 & Zi1[, 2]==0)
    
    #0～1階層目の対数尤度
    LL0 <- sum(log(phi0[wd[Zi1[, 1]==0]]))
    LL1 <- sum(log((t(phi1)[wd[index], ] * Zi21[index, ]) %*% rep(1, k1)))
    
    #2～3階層目の対数尤度
    LL2 <- c(); LL3 <- c()
    for(i in 1:k1){
      index2 <- which(Zi21[, i]==1 & (Zi_path*(1-Zi1[, 3]))==1)
      LLi2 <- as.numeric((t(phi2[[i]])[wd[index2], ] * Zi22_list[[i]][index2, ]) %*% rep(1, k2[i]))
      
      LL2 <- c(LL2, sum(log(LLi2[LLi2!=0])))
      for(j in 1:k2[i]){
        index3 <- which(Zi22_list[[i]][, j]==1 & Zi1[, 3]==1)
        LLi3 <- (t(phi3[[i]][[j]])[wd[index3], , drop=FALSE] * (Zi23_list[[i]][[j]])[index3, ]) %*% rep(1, ncol(Zi23_list[[i]][[j]]))
        LL3 <- c(LL3, sum(log(LLi3[LLi3!=0])))
      }
    }
    LL <- LL0 + LL1 + sum(LL2) + sum(LL3)   #対数尤度の総和
    
    #サンプリング結果の表示
    print(rp)
    print(c(LL, LLst))
    print(c(LL0, LL1, sum(LL2), sum(LL3)))
    print(round(rbind(theta1, thetat1), 3))
    print(round(rbind(print(rep(1:k1, k2)), unlist(gamma3), unlist(gammat3)), 3))
  }
}

