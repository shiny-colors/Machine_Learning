#####Hidden Marcov language model#####
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

#set.seed(92483)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k1 <- 5   #��ʌ�̕��@������
k2 <- 7   #�@�\��̕i��������
d <- 2000   #������
v1 <- 750   #��ʌꐔ
v2 <- 50   #�@�\�ꐔ
s <- rpois(d, 13)   #���͐�
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
a <- sum(s)   #�����͐�
w <- rpois(a, 16)   #���͂�����̒P�ꐔ
w[w < 5] <- ceiling(runif(sum(w < 5), 5, 10))
f <- sum(w)   #���P�ꐔ

#����ID�̐ݒ�
u_id <- rep(1:d, s)
t_id <- c()
for(i in 1:d){t_id <- c(t_id, 1:s[i])}
words <- as.numeric(tapply(w, u_id, sum))

#���͋�؂�̃x�N�g�����쐬
ID_d <- rep(1:d, words)
td_d <- c()
for(i in 1:a){
  td_d <- c(td_d, 1:w[i])
}
nd_d <- rep(1:a, w)
x_vec <- rep(0, f)
x_vec[c(1, cumsum(w[-a])+1)] <- 1

#�C���f�b�N�X��ݒ�
s_list <- list()
vec_list <- list()
for(i in 1:a){
  if(i%%1000==0){
    print(i)
  }
  s_list[[i]] <- which(nd_d==i)
  vec_list[[i]] <- rep(1, length(s_list[[i]]))
}


##�p�����[�^�̐ݒ�
#�}���R�t���ڍs��̃p�����[�^
alpha01 <- rep(2.0, k1+1)
alpha02 <- rep(0.2, k1+1)
alpha03 <- rep(0.2, k2)
theta1 <- thetat1 <- as.numeric(extraDistr::rdirichlet(1, alpha01))
theta2 <- thetat2 <- extraDistr::rdirichlet(k2+1, alpha02)
theta3 <- thetat3 <- extraDistr::rdirichlet(k1, alpha03)

#�P��o�����̃p�����[�^
alpha11 <- rep(0.1, v1)
alpha12 <- rep(0.1, v2)
phi1 <- phit1 <- extraDistr::rdirichlet(k1+1, alpha11)
phi2 <- phit2 <- extraDistr::rdirichlet(k2, alpha12)


##���ݐ��ڂƊϑ��f�[�^�̐���
Z1_list <- list() 
Z2_list <- list()
WX1 <- matrix(0, nrow=a, ncol=v1)
WX2 <- matrix(0, nrow=a, ncol=v2)
word_list <- list()

for(i in 1:a){
  if(i%%1000==0){
    print(i)
  }
  
  #�f�[�^�̐ݒ�
  n <- w[i]
  z1 <- matrix(0, nrow=n, ncol=k1+1)
  z2 <- matrix(0, nrow=n, ncol=k2)
  z1_vec <- rep(0, n)
  z2_vec <- rep(0, n)
  words1 <- matrix(0, nrow=n, ncol=v1)
  words2 <- matrix(0, nrow=n, ncol=v2)
  
  #�P�ꂲ�Ƃɕ��@����я����𐶐�
  for(j in 1:n){
    
    #���͂�1�P��ڂ̕��@�𐶐�
    if(j==1){
      z1[j, ] <- rmnom(1, 1, theta1)
      z1_vec[j] <- as.numeric(z1[j, ] %*% 1:(k1+1))
      words1[j, ] <- rmnom(1, 1, phi1[z1_vec[j], ])
    }
    
    #2�P��ڈȍ~�̕i���y�я����𐶐�
    if(j >= 2){
      
      #���@�ƒP��𐶐�
      if(z1_vec[j-1]==(k1+1)){
        z1[j, ] <- rmnom(1, 1, theta2[k2+1, ])
        z1_vec[j] <- as.numeric(z1[j, ] %*% 1:(k1+1))
        words1[j, ] <- rmnom(1, 1, phi1[z1_vec[j], ])
        
      } else if(z2_vec[j-1]!=0){
        z1[j, ] <- rmnom(1, 1, theta2[z2_vec[j-1], ])
        z1_vec[j] <- as.numeric(z1[j, ] %*% 1:(k1+1))
        words1[j, ] <- rmnom(1, 1, phi1[z1_vec[j], ])
    
      #�����ƒP��𐶐�  
      } else {
        z2[j, ] <- rmnom(1, 1, theta3[z1_vec[j-1], ])
        z2_vec[j] <- as.numeric(z2[j, ] %*% 1:k2)
        words2[j, ] <- rmnom(1, 1, phi2[z2_vec[j], ])
      }
    }
  }
  #�f�[�^���i�[
  Z1_list[[i]] <- z1_vec
  Z2_list[[i]] <- z2_vec
  WX1[i, ] <- colSums(words1)
  WX2[i, ] <- colSums(words2)
  word_list[[i]] <- as.numeric(words1 %*% 1:v1 + words2 %*% 1:v2)
}



##���X�g��ϊ�
Z1 <- unlist(Z1_list)
Z2 <- unlist(Z2_list)
word_vec02 <- word_vec01 <- unlist(word_list)
word_vec1 <- word_vec01 * as.numeric(Z1 > 0)
word_vec2 <- word_vec02 * as.numeric(Z2 > 0)
storage.mode(WX1) <- "integer"
storage.mode(WX2) <- "integer"
sparse_data1 <- as(WX1, "CsparseMatrix")
sparse_data2 <- as(WX2, "CsparseMatrix")

#�X�p�[�X�s����쐬
index_part <- which(word_vec2 > 0)   #�����̃C���f�b�N�X
i1 <- which(word_vec1 > 0); j1 <- word_vec1[word_vec1 > 0]; x1 <- rep(1, length(word_vec1[word_vec1 > 0]))
i2 <- which(word_vec2 > 0); j2 <- word_vec2[word_vec2 > 0]; x2 <- rep(1, length(word_vec2[word_vec2 > 0]))
sparse_word1 <- sparseMatrix(i1, j1, x=x1, dims=c(f, v1))
sparse_word2 <- sparseMatrix(i2, j2, x=x2, dims=c(f, v2))
word_data1 <- as.matrix(sparse_word1)
word_data2 <- as.matrix(sparse_word2)
storage.mode(word_data1) <- "integer"
storage.mode(word_data2) <- "integer"


####�}���R�t�A�������e�J�����@��Hidden Marcov language model�𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##���O���z�̐ݒ�
alpha01 <- 1
alpha02 <- 1
beta01 <- 0.1
beta02 <- 0.1

##�p�����[�^�̐^�l
phi1 <- phit1
phi2 <- phit2
theta1 <- thetat1
theta2 <- thetat2
theta3 <- thetat3

##�����l�̐ݒ�
#�}���R�t���ڍs��̏����l
theta1 <- as.numeric(extraDistr::rdirichlet(1, rep(1, k1+1)))
theta2 <- extraDistr::rdirichlet(k2+1, rep(1.0, k1+1))
theta3 <- extraDistr::rdirichlet(k1, rep(1.0, k2))

#�P��o�����̏����l
phi1 <- extraDistr::rdirichlet(k1+1, rep(0.5, v1))
phi2 <- extraDistr::rdirichlet(k2, rep(0.5, v2))


##�p�����[�^�̊i�[�p�z��
THETA1 <- matrix(0, nrow=R/keep, ncol=k1+1)
THETA2 <- array(0, dim=c(k2+1, k1+1, R/keep))
THETA3 <- array(0, dim=c(k2, k1, R/keep))
PHI1 <- array(0, dim=c(k1+1, v1, R/keep))
PHI2 <- array(0, dim=c(k2, v2, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k1+1)
SEG2 <- matrix(0, nrow=f, ncol=k2)

##�f�[�^�̐ݒ�
#��ʌꂪ�A�����Ă��镔�����C���f�b�N�X��
index1 <- which(as.numeric(word_vec1 > 0)==1)[-1]
index2 <- index1+1
index_phrase <- subset(index1, word_vec1[index2] > 0 & x_vec[index2]==0)

#�@�\��ƈ�ʌ�̃C���f�b�N�X���쐬
index_general <- which(rowSums(sparse_word1) > 0)
index_function <- which(rowSums(sparse_word2) > 0)

#�P��̃C���f�b�N�X���쐬
'%!in%' <- function(x,y)!('%in%'(x,y))
max_pt <- max(w)
index_t11 <- which(td_d==1)
index_t21 <- list()
index_t22 <- list()
index_t21_g1 <- list()
index_t22_g1 <- list()
index_t21_g2 <- list()
index_t22_g2 <- list()
index_t21_f <- list()
index_t22_f <- list()

for(j in 2:max_pt){
  index_t21[[j]] <- which(td_d==j)-1
  index_t22[[j]] <- index_t21[[j]]+1
  index_t21[[j]][rowSums(sparse_word1[index_t22[[j]], , drop=FALSE]) > 0]
  index_t21_g1[[j]] <- index_t21[[j]][rowSums(sparse_word1[index_t22[[j]], , drop=FALSE]) > 0]
  index_t22_g1[[j]] <- index_t22[[j]][rowSums(sparse_word1[index_t22[[j]], , drop=FALSE]) > 0]
  index_g1 <- index_t21[[j]][rowSums(sparse_word1[index_t22[[j]], , drop=FALSE]) > 0] 
  index_t21_g1[[j]] <- subset(index_g1, index_g1 %!in% index_phrase) 
  index_t22_g1[[j]] <- index_t21_g1[[j]]+1 
  index_t21_g2[[j]] <- subset(index_g1, index_g1 %in% index_phrase)
  index_t22_g2[[j]] <- index_t21_g2[[j]]+1 
  index_t22_f[[j]] <- index_t22[[j]][rowSums(sparse_word2[index_t22[[j]], , drop=FALSE]) > 0]
  index_t21_f[[j]] <- index_t22_f[[j]]-1
}

##�ΐ��ޓx�̊�l
LLst1 <- sum(sparse_word1[index_general, ] %*% t(log(extraDistr::rdirichlet(1, colSums(sparse_word1) + beta01))))
LLst2 <- sum(sparse_word2[index_function, ] %*% t(log(extraDistr::rdirichlet(1, colSums(sparse_word2) + beta02))))
LLst <- LLst1 + LLst2


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  #�f�[�^�̊i�[�p�z��
  Zi1 <- matrix(0, nrow=f, ncol=k1+1)
  z1_vec <- rep(0, f)
  Zi2 <- matrix(0, nrow=f, ncol=k2)
  z2_vec <- rep(0, f)
  rf01 <- matrix(0, nrow=k2+1, ncol=k1+1)
  rf02 <- matrix(0, nrow=k1, ncol=k2)
  
  ##1�P��ڂ̐��ݕϐ��𐶐�
  #���ݕϐ��̃p�����[�^��ݒ�
  n <- length(index_t11)
  z_par1 <- sparse_word1[index_t11, ] %*% t(phi1[-(k1+1), ]) * matrix(theta1[-(k1+1)], nrow=n, ncol=k1, byrow=T)
  z_rate1 <- z_par1 / rowSums(z_par1)

  #�������z������ݕϐ��𐶐�
  Zi1[index_t11, -(k1+1)] <- rmnom(n, 1, z_rate1)
  Zi1[index_t21_g2[[2]], ] <- 0; Zi1[index_t21_g2[[2]], k1+1] <- 1   #��ʌ�̘A���������m��
  z1_vec[index_t11] <- Zi1[index_t11, ] %*% 1:(k1+1)
  
  
  ##2�P��ڈȍ~�̐��ݕϐ��𒀎�����
  for(j in 2:max_pt){
    
    ##�@�\��̐��ݕϐ��𐶐�
    if(length(index_t21_f[[j]]) > 0){
      
      #���ݕϐ��̃p�����[�^��ݒ�
      index1 <- index_t21_f[[j]]
      index2 <- index_t22_f[[j]]
      z_par2 <- sparse_word2[index2, ] %*% t(phi2) * theta3[z1_vec[index1], ]   #���ݕϐ��̃p�����[�^
      z_rate2 <- z_par2 / rowSums(z_par2)
      
      #�������z������ݕϐ��𐶐�
      Zi2[index2, ] <- rmnom(length(index1), 1, z_rate2)
      z2_vec[index2] <- as.numeric(Zi2[index2, ] %*% 1:k2)
      
      #�}���R�t���ڍs����X�V
      rf02 <- rf02 + t(Zi1[index_t21[[j]], -(k1+1), drop=FALSE]) %*% Zi2[index_t22[[j]], , drop=FALSE]   
    }
    
    ##��t���[�Y�̈�ʌ�̐��ݕϐ��𐶐�
    if(length(index_t21_g1[[j]]) > 0){
      
      #���ݕϐ��̃p�����[�^��ݒ�
      index3 <- index_t21_g1[[j]]
      index4 <- index_t22_g1[[j]]
      n <- length(index4)
      z_par1 <- sparse_word1[index4, , drop=FALSE] %*% t(phi1[-(k1+1), ]) * theta2[z2_vec[index3], -(k1+1)]
      z_rate1 <- z_par1 / rowSums(z_par1)
      
      #�������z������ݕϐ��𐶐�
      Zi1[index4, -(k1+1)] <- rmnom(n, 1, z_rate1)
      if(j < max_pt){
        Zi1[index_t21_g2[[j+1]], ] <- 0; Zi1[index_t21_g2[[j+1]], k1+1] <- 1   #��ʌ�̘A���������m��
      }
      z1_vec[index4] <- Zi1[index4, ] %*% 1:(k1+1)
    }
    
    ##�t���[�Y�̈�ʌ�̐��ݕϐ��𐶐�
    if(length(index_t21_g2[[j]]) > 0){
      
      #���ݕϐ��̃p�����[�^��ݒ�
      index5 <- index_t21_g2[[j]]
      index6 <- index_t22_g2[[j]]
      n <- length(index5)
      z_par1 <- sparse_word1[index6, ] %*% t(phi1[-(k1+1), ]) * matrix(theta2[k2+1, (-k1+1)], nrow=n, ncol=k1, byrow=T)
      z_rate1 <- z_par1 / rowSums(z_par1)
      
      #�������z������ݕϐ��𐶐�
      Zi1[index6, -(k1+1)] <- rmnom(n, 1, z_rate1)
      if(j < max_pt){
        Zi1[index_t21_g2[[j+1]], ] <- 0; Zi1[index_t21_g2[[j+1]], k1+1] <- 1   #��ʌ�̘A���������m��
      }
      z1_vec[index6] <- Zi1[index6, ] %*% 1:(k1+1)
    }
    #�}���R�t���ڍs����X�V
    rf01 <- rf01 + t(cbind(Zi2[index_t21[[j]], , drop=FALSE], Zi1[index_t21[[j]], k1+1])) %*% Zi1[index_t22[[j]], , drop=FALSE]
  }
  
  
  ##�}���R�t���ڍs��̃p�����[�^���X�V
  #�f�B���N�����z�̃p�����[�^��ݒ�
  par1 <- colSums(Zi1[index_t11, ]) + alpha01
  rf1 <- rf01 + alpha01
  rf2 <- rf02 + alpha02
  
  #�f�B���N�����z����P��o�������T���v�����O
  theta1 <- as.numeric(extraDistr::rdirichlet(1, par1))
  theta2 <- extraDistr::rdirichlet(k2+1, rf1)
  theta3 <- extraDistr::rdirichlet(k1, rf2)
  
  
  ##�P��o���m���̃p�����[�^���X�V
  #�f�B���N�����z�̃p�����[�^��ݒ�
  wf1 <- t(as.matrix(t(sparse_word1[index_general, ]) %*% Zi1[index_general, ])) + beta01
  wf2 <- t(as.matrix(t(sparse_word2[index_function, ]) %*% Zi2[index_function, ])) + beta02
  phi1 <- extraDistr::rdirichlet(k1+1, wf1)
  phi2 <- extraDistr::rdirichlet(k2, wf2)
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    mkeep <- 1
    THETA1[mkeep, ] <- theta1
    THETA2[, , mkeep] <- theta2
    THETA3[, , mkeep] <- theta3
    PHI1[, , mkeep] <- phi1
    PHI2[, , mkeep] <- phi2
    
    #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
    if(rp%%keep==0 & rp >= burnin){
      SEG1 <- SEG1 + Zi1
      SEG2 <- SEG2 + Zi2
    }
  
    #�T���v�����O���ʂ��m�F
    if(rp%%disp==0){
      #�ΐ��ޓx���v�Z
      LL1 <- sum(sparse_word1[index_general, ] %*% t(log(phi1)) * Zi1[index_general, ])
      LL2 <- sum(sparse_word2[index_function, ] %*% t(log(phi2)) * Zi2[index_function, ])
      
      #�T���v�����O���ʂ̕\��
      print(rp)
      print(c(LL1+LL2, LLst))
      print(round(rbind(theta1, thetat1), 3))
      print(round(cbind(theta3, thetat3), 3))
      print(round(cbind(phi1[, 1:10], phit1[, 1:10]), 3))
      print(round(cbind(phi2[, 1:10], phit2[, 1:10]), 3))
    }
  }
}
