#####�x�C�W�A���m���I�����u���b�N���f��#####
library(MASS)
library(Matrix)
library(bayesm)
library(matrixStats)
library(MCMCpack)
library(gtools)
library(extraDistr)
library(reshape2)
library(qrmtools)
library(slfm)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
N <- 10000   #���[�U�[��
K <- 3000   #�A�C�e����
seg_u <- 8   #���[�U�[�̃Z�O�����g��
seg_i <- 7   #�A�C�e���̃Z�O�����g��

##�p�����[�^�ƃZ�O�����g�𔭐�������
#���[�U�[�Z�O�����g�𔭐�
alpha01 <- rep(25, seg_u)
pi01 <- extraDistr::rdirichlet(1, alpha01)
z01 <- t(rmultinom(N, 1, pi01))
z01_vec <- as.numeric(z01 %*% 1:seg_u)
z1_cnt <- as.numeric(table(z01_vec))
mix1 <- colMeans(z01)   #�������̐^�l


#�A�C�e���Z�O�����g�̔���
alpha02 <- rep(25, seg_i)
pi02 <- extraDistr::rdirichlet(1, alpha02)
z02 <- t(rmultinom(K, 1, pi02))
z02_vec <- as.numeric(z02 %*% 1:seg_i)
z2_cnt <- as.numeric(table(z02_vec))
mix2 <- colMeans(z02)   #�������̐^�l


#�ϑ��ϐ��̃p�����[�^�̐ݒ�
#���[�U�[�Z�O�����g�~�A�C�e���Z�O�����g�̃x�[�^���O���z�̃p�����[�^�𔭐�
hist(rbeta(10000, 1.0, 6.5), col="grey", breaks=25, main="�x�[�^���z����̗���", xlab="�p�����[�^")
theta0 <- matrix(rbeta(seg_u*seg_i, 1.0, 6.5), nrow=seg_u, ncol=seg_i)
round(theta0, 3)


##�x���k�[�C���z����ϑ��s��𔭐�������
#�x���k�[�C���z���狤�N�s��𐶐�
Data <- matrix(0, nrow=N, ncol=K)

for(i in 1:seg_u){
  print(i)
  for(j in 1:seg_i){
    n <- z1_cnt[i] * z2_cnt[j]
    Data[z01_vec==i, z02_vec==j] <- matrix(rbinom(n, 1, theta0[i, j]), nrow=z1_cnt[i], ncol=z2_cnt[j])
  }
}
Data_T <- t(Data)
storage.mode(Data) <- "integer"
storage.mode(Data_T) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")
sparse_data_T <- as(Data_T, "CsparseMatrix")  
gc(); gc()

#�������z�̃p�����[�^��ݒ�
thetat <- matrix(0, nrow=seg_u, ncol=seg_i)
for(i in 1:seg_u){
  n <- sum(sparse_data[z01_vec==i, ])
  for(j in 1:seg_i){
    thetat[i, j] <- sum(sparse_data[z01_vec==i, z02_vec==j]) / n
  }
}

phit <- matrix(0, nrow=seg_i, nco=seg_u)
for(i in 1:seg_i){
  n <- sum(sparse_data[, z02_vec==i])
  for(j in 1:seg_u){
    phit[i, j] <- sum(sparse_data[z01_vec==j, z02_vec==i]) / n
  }
}


####�}���R�t�A�������e�J�����@�Ŋm���I�u���b�N���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##�萔�̐ݒ�
const_user <- lfactorial(n_user)
const_item <- lfactorial(n_item)

##���O���z�̐ݒ�
alpha1 <- matrix(5.0, nrow=seg_u, ncol=seg_i)   #���[�U�[�Z�O�����g�̃f�B�N�������O���z
alpha2 <- matrix(5.0, nrow=seg_i, ncol=seg_u)   #�A�C�e���Z�O�����g�̃f�B�N�������O���z

##�����l�̐ݒ�
#�������̏����l
r1 <- rep(1/seg_u, seg_u)
r2 <- rep(1/seg_i, seg_i)

#�u���b�N���Ƃ̃p�����[�^�̏����l
index_u <- floor(seq(1, N, length=seg_u+1))
index_i <- floor(seq(1, K, length=seg_i+1))
sortlist1 <- order(rowSums(sparse_data))
sortlist2 <- order(colSums(sparse_data), decreasing=TRUE)


#�N���X�����̏����l��ݒ�
z1 <- rep(0, N)
z2 <- rep(0, K)

for(i in 1:(length(index_u)-1)){
  #���[�U�[�̃N���X�����̃C���f�b�N�X��ݒ�
  index1 <- sortlist1[index_u[i]:index_u[i+1]]
  z1[index1] <- i
  
  for(j in 1:(length(index_i)-1)){
    #�A�C�e���̃N���X�����̃C���f�b�N�X��ݒ�
    index2 <- sortlist2[index_i[j]:index_i[j+1]]
    z2[index2] <- j
  }
}

#�Z�O�����g�C���f�b�N�X���쐬
index1 <- list()
index2 <- list()
item_vec <- matrix(0, nrow=K, ncol=seg_i)
user_vec <- matrix(0, nrow=N, ncol=seg_u)
n1 <- c()
n2 <- c()

for(i in 1:seg_u){
  index1[[i]] <- which(z1==i)
  user_vec[index1[[i]], i] <- 1
  n1 <- c(n1, sum(sparse_data[index1[[i]], ]))
}
for(j in 1:seg_i){
  index2[[j]] <- which(z2==j)
  item_vec[index2[[j]], j] <- 1
  n2 <- c(n2, sum(sparse_data[, index2[[j]]]))
}

#�p�����[�^�̏����l��ݒ�
#�A�C�e���Z�O�����g�̏����p�����[�^
oldtheta <- matrix(0, nrow=seg_u, ncol=seg_i) 
for(i in 1:seg_u){
  for(j in 1:(seg_i-1)){
    freq <- sum(sparse_data[index1[[i]], index2[[j]]])
    oldtheta[i, j] <- freq / n1[i]
  }
}
oldtheta[, seg_i] <- 1 - rowSums(oldtheta)

#�A�C�e���Z�O�����g�̏����p�����[�^
oldphi <- matrix(0, nrow=seg_i, ncol=seg_u)
for(i in 1:seg_i){
  for(j in 1:(seg_u-1)){
    freq <- sum(sparse_data[index1[[j]], index2[[i]]])
    oldphi[i, j] <- freq / n2[i]
  }
}
oldphi[, seg_u] <- 1 - rowSums(oldphi)


#���[�U�[�A�A�C�e�����Ƃ̍w����
n_user <- rowSums(sparse_data)
n_item <- colSums(sparse_data)

##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(seg_u, seg_i, R/keep))
PHI <- array(0, dim=c(seg_i, seg_u, R/keep))
SEG1 <- matrix(0, nrow=R/keep, ncol=N)
SEG2 <- matrix(0, nrow=R/keep, ncol=K)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()



####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���[�U�[�̃Z�O�����g�����𐶐�
  y1 <- as.matrix(sparse_data %*% item_vec)   #�A�C�e���Z�O�����g���Ƃ̍w���p�x
  LLi1 <- y1 %*% log(t(oldtheta))   #�Z�O�����g���Ƃɑ������z�̖ޓx���v�Z
  
  #logsumexp�̖ޓx���v�Z
  LLi_max <- rowMaxs(LLi1)
  r_matrix <- matrix(r1, nrow=N, ncol=seg_u, byrow=T)
  
  #�����m���̃p�����[�^��ݒ�
  expl <- r_matrix * exp(LLi1 - LLi_max)
  z1_rate <- expl / rowSums(expl)   #�Z�O�����g�����m��

  
  #�������z����Z�O�����g�����𐶐�
  Z1 <- rmnom(N, 1, z1_rate)
  z1 <- as.numeric(Z1 %*% 1:seg_u)
  
  #���������X�V
  z_sums <- colSums(Z1) + 1
  r1 <- z_sums / sum(z_sums)
  
  
  ##�A�C�e���̃Z�O�����g�����𐶐�
  y2 <- as.matrix(sparse_data_T %*% user_vec)   #���[�U�Z�O�����g���Ƃ̍w���p�x
  LLi2 <- y2 %*% t(log(oldphi))   #�Z�O�����g���Ƃɑ������z�̖ޓx���v�Z

  #logsumexp�̖ޓx���v�Z
  LLi_max <- rowMaxs(LLi2)
  r_matrix <- matrix(r2, nrow=K, ncol=seg_i, byrow=T)
  
  #�����m���̃p�����[�^��ݒ�
  expl <- r_matrix * exp(LLi2 - LLi_max)
  z2_rate <- expl / rowSums(expl)   #�Z�O�����g�����m��
  
  #�������z����Z�O�����g�����𐶐�
  Z2 <- rmnom(K, 1, z2_rate)
  z2 <- as.numeric(Z2 %*% 1:seg_i)
  
  #���������X�V
  z_sums <- colSums(Z2) + 1
  r2 <- z_sums / sum(z_sums)
  
  
  ##�Z�O�����g�C���f�b�N�X���쐬
  index1 <- list()
  index2 <- list()
  item_vec <- matrix(0, nrow=K, ncol=seg_i)
  user_vec <- matrix(0, nrow=N, ncol=seg_u)
  
  for(i in 1:seg_u){
    index1[[i]] <- which(z1==i)
    user_vec[index1[[i]], i] <- 1
  }
  for(j in 1:seg_i){
    index2[[j]] <- which(z2==j)
    item_vec[index2[[j]], j] <- 1
  }
  
  ##���[�U�[����уA�C�e���̃p�����[�^���T���v�����O
  #�f�B�N�������z���烆�[�U�[�Z�O�����g�̃p�����[�^���T���v�����O
  freq_user <- matrix(0, nrow=seg_u, ncol=seg_i)
  for(i in 1:seg_u){
    for(j in 1:seg_i){
      freq_user[i, j] <- sum(sparse_data[index1[[i]], index2[[j]]])
    }
  }
  oldtheta <- extraDistr::rdirichlet(seg_u, freq_user + alpha1)
  
  #�f�B�N�������z����A�C�e���Z�O�����g�̃p�����[�^���T���v�����O
  freq_item <- matrix(0, nrow=seg_i, ncol=seg_u)
  for(i in 1:seg_i){
    for(j in 1:seg_u){
      freq_item[i, j] <- sum(sparse_data[index1[[j]], index2[[i]]])
    }
  }
  oldphi <- extraDistr::rdirichlet(seg_i, freq_item + alpha2)
  
  ##�T���v�����O���ʂ̊i�[�ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta 
    PHI[, , mkeep] <- oldphi
    SEG1[mkeep, ] <- z1 
    SEG2[mkeep, ] <- z2 
    
    if(rp%%disp==0){
      print(rp)
      print(round(rbind(r1, mix1), 3))
      print(round(rbind(r2, mix2), 3))
      print(round(cbind(oldtheta, thetat), 3))
      print(round(cbind(oldphi, phit), 3))
    }
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 1000/keep   #�o�[���C������ 
r <- R/keep   #�T���v�����O�̍ŏI�s

##�T���v�����O���ʂ̉���
matplot(t(THETA[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(THETA[2, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(THETA[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(THETA[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(THETA[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(PHI[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(PHI[2, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(PHI[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")
matplot(t(PHI[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^", main="�������z�̃p�����[�^")

##�p�����[�^�̐���l
#���[�U�[�Z�O�����g�̐���l�̗v�񓝌v��
round(cbind(apply(THETA[, , burnin:r], c(1, 2), mean), thetat), 3)   #���㕽��
round(apply(THETA[, , burnin:r], c(1, 2), function(x) quantile(x, 0.025)), 3)   #����5�����ʓ_
round(apply(THETA[, , burnin:r], c(1, 2), function(x) quantile(x, 0.975)), 3)   #����95�����ʓ_
round(apply(THETA[, , burnin:r], c(1, 2), sd), 3)   #����W���΍�

#�A�C�e���Z�O�����g�̐���l�̗v�񓝌v��
round(cbind(apply(PHI[, , burnin:r], c(1, 2), mean), phit), 3)   #���㕽��
round(apply(PHI[, , burnin:r], c(1, 2), function(x) quantile(x, 0.025)), 3)   #����5�����ʓ_
round(apply(PHI[, , burnin:r], c(1, 2), function(x) quantile(x, 0.975)), 3)   #����95�����ʓ_
round(apply(PHI[, , burnin:r], c(1, 2), sd), 3)   #����W���΍�

##�Z�O�����g�����Ƌ��N�֌W������
#���[�U�[�Z�O�����g�̊������m��
seg1 <- rep(0, N)
for(i in 1:N){
  x <- table(SEG1[burnin:r, i])
  seg1[i] <- as.numeric(names(which.max(x)))
}

#�A�C�e���Z�O�����g�̊������m��
seg2 <- rep(0, K)
for(i in 1:K){
  x <- table(SEG2[burnin:r, i])
  seg2[i] <- as.numeric(names(which.max(x)))
}

#�f�[�^���Z�O�����g���ʕ��ёւ�
index_seg1 <- order(seg1)
index_seg2 <- order(seg2)
Block_Data <- Data[index_seg1, index_seg2]   #���f�[�^���Z�O�����g���ɕ��ёւ���

#���ёւ����u���b�N�f�[�^������
plot_matrix(Data, standardize.rows=FALSE, reorder.rows=FALSE, reorder.cols=FALSE, high.contrast=TRUE)
plot_matrix(Block_Data, standardize.rows=FALSE, reorder.rows=FALSE, reorder.cols=FALSE, high.contrast=TRUE)

