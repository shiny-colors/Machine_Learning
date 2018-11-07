#####Hidden Marcov Stochastic Block Multinomial model#####
library(MASS)
library(Matrix)
library(bayesm)
library(MCMCpack)
library(gtools)
library(extraDistr)
library(matrixStats)
library(reshape2)
library(qrmtools)
library(slfm)
library(HMM)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

#set.seed(318)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
hh <- 3000   #���[�U�[��
pt <- rpois(hh, 20)   #�ϑ�����
max_pt <- max(pt)
hhpt <- sum(pt)   #�����R�[�h��
item <- 2000   #�A�C�e����
seg_u <- 8   #���[�U�[�̃Z�O�����g��
seg_i <- 7   #�A�C�e���̃Z�O�����g��

#ID��ݒ�
u_id <- rep(1:hh, pt)
t_id <- c()
for(i in 1:hh){
  t_id <- c(t_id, 1:pt[i])
}

#�C���f�b�N�X���쐬
u_index <- list()
for(i in 1:hh){u_index[[i]] <- which(u_id==i)}

##�p�����[�^�ƃZ�O�����g�𐶐�
#�p�����[�^�𐶐�
alpha01 <- rep(10, seg_u)
alpha02 <- matrix(2.0, nrow=seg_u, ncol=seg_u)
diag(alpha02) <- runif(seg_u, 10, 16)
gamma01 <- as.numeric(extraDistr::rdirichlet(1, alpha01))
gamma02 <- extraDistr::rdirichlet(seg_u, alpha02)

##���[�U�[�Z�O�����g�𐶐�
z01 <- matrix(0, nrow=hhpt, ncol=seg_u)
z01_vec <- rep(0, hhpt)
for(i in 1:hh){
  for(j in 1:pt[i]){
    if(j==1){
      #1���ڂ̃Z�O�����g�𐶐�
      z01[u_index[[i]][j], ] <- rmnom(1, 1, gamma01)
      z01_vec[u_index[[i]][j]] <- as.numeric(z01[u_index[[i]][j], ] %*% 1:seg_u)
      
    } else {
      
      #2���ڈȍ~�̃Z�O�����g�𐶐�
      z01[u_index[[i]][j], ] <- rmnom(1, 1, gamma02[z01_vec[u_index[[i]][j-1]], ])
      z01_vec[u_index[[i]][j]] <- as.numeric(z01[u_index[[i]][j], ] %*% 1:seg_u)
    }
  }
}
z1_cnt <- colSums(z01)

##�A�C�e���Z�O�����g�𐶐�
alpha03 <- rep(20, seg_i)
omega01 <- as.numeric(extraDistr::rdirichlet(1, alpha03))
z02 <- rmnom(item, 1, omega01)
z02_vec <- as.numeric(z02 %*% 1:seg_i)
z2_cnt <- colSums(z02)


##�ϑ����f���̃p�����[�^�𐶐�
#���[�U�[�Z�O�����g�~�A�C�e���Z�O�����g�̃x�[�^���O���z�̃p�����[�^�𔭐�
hist(rbeta(10000, 0.4, 4.5), col="grey", breaks=25, main="�x�[�^���z����̗���", xlab="�p�����[�^")
theta0 <- matrix(rbeta(seg_u*seg_i, 0.4, 4.5), nrow=seg_u, ncol=seg_i)
round(theta0, 3)

##�x���k�[�C���z����ϑ��s��𔭐�������
#�x���k�[�C���z���狤�N�s��𐶐�
Data <- matrix(0, nrow=hhpt, ncol=item)

for(i in 1:seg_u){
  print(i)
  for(j in 1:seg_i){
    n <- z1_cnt[i] * z2_cnt[j]
    z1_cnt[i]
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


####�}���R�t�A�������e�J�����@��HM Sparse Stochastic Block model�𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000 
keep <- 2
disp <- 10
iter <- 0
sbeta <- 1.5

##���O���z�̐ݒ�
tau1 <- rep(1, seg_u)
tau2 <- matrix(1, nrow=seg_u, ncol=seg_u)
diag(tau2) <- 5
tau3 <- rep(1, seg_i)
alpha1 <- matrix(1, nrow=seg_u, ncol=seg_i)   #���[�U�[�Z�O�����g�̃f�B�N�������O���z
alpha2 <- matrix(1, nrow=seg_i, ncol=seg_u)   #�A�C�e���Z�O�����g�̃f�B�N�������O���z


##�����l�̐ݒ�
#�������̏����l
alpha01 <- rep(20, seg_u)
alpha02 <- matrix(4, nrow=seg_u, ncol=seg_u)
alpha03 <- rep(20, seg_i)
diag(alpha02) <- 20
gamma1 <- extraDistr::rdirichlet(1, alpha01)
gamma2 <- extraDistr::rdirichlet(seg_u, alpha02)
omega1 <- extraDistr::rdirichlet(1, alpha03)

#�u���b�N���Ƃ̃p�����[�^�̏����l
index_u <- floor(seq(1, hhpt, length=seg_u+1))
index_i <- floor(seq(1, item, length=seg_i+1))
sortlist1 <- order(rowSums(sparse_data))
sortlist2 <- order(colSums(sparse_data), decreasing=TRUE)

#�N���X�����̏����l��ݒ�
z1 <- rep(0, hhpt)
z2 <- rep(0, item)

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
user_vec <- matrix(0, nrow=hhpt, ncol=seg_u)
item_vec <- matrix(0, nrow=item, ncol=seg_i)
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
oldtheta <- (oldtheta + 0.00001) / rowSums(oldtheta + 0.00001)

#�A�C�e���Z�O�����g�̏����p�����[�^
oldphi <- matrix(0, nrow=seg_i, ncol=seg_u)
for(i in 1:seg_i){
  for(j in 1:(seg_u-1)){
    freq <- sum(sparse_data[index1[[j]], index2[[i]]])
    oldphi[i, j] <- freq / n2[i]
  }
}
oldphi[, seg_u] <- 1 - rowSums(oldphi)
oldphi <- (oldphi + 0.00001) / rowSums(oldphi + 0.00001)

#���[�U�[�A�A�C�e�����Ƃ̍w����
n_user <- rowSums(sparse_data)
n_item <- colSums(sparse_data)

##�p�����[�^�̊i�[�p�z��
THETA <- array(0, dim=c(seg_u, seg_i, R/keep))
PHI <- array(0, dim=c(seg_i, seg_u, R/keep))
GAMMA1 <- matrix(0, nrow=R/keep, ncol=seg_u)
GAMMA2 <- array(0, dim=c(seg_u, seg_u, R/keep))
OMEGA1 <- matrix(0, nrow=R/keep, ncol=seg_i)
SEG1 <- matrix(0, hhpt, ncol=seg_u)
SEG2 <- matrix(0, item, ncol=seg_i)
storage.mode(SEG1) <- "integer"
storage.mode(SEG2) <- "integer"
gc(); gc()


##�C���f�b�N�X���쐬
max_pt <- max(pt)
index_t11 <- which(t_id==1)
index_t21 <- list()
index_t22 <- list()
n_time <- rep(0, max_pt)
n_time[1] <- length(index_t11)
for(j in 2:max_pt){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
  n_time[j] <- length(index_t22[[j]])
}


####�}���R�t�A�������e�J�����@�Ńp�����[�^���T���v�����O####
for(rp in 1:R){

  ##���[�U�[���Ƃ̃Z�O�����g�����𐶐�
  #���[�U�[�̃Z�O�����g�������Ƃ̖ޓx�𐄒�
  y1 <- as.matrix(sparse_data %*% item_vec)   #�A�C�e���Z�O�����g���Ƃ̍w���p�x
  LLi0 <- y1 %*% t(log(oldtheta))   #�Z�O�����g���Ƃ̑������z�̑ΐ��ޓx
  LLi_max <- rowMaxs(LLi0)
  LLi1 <- exp(LLi0 - LLi_max)   #�ޓx�ɕϊ�
  
  #�Z�O�����g�����m���̐���ƃZ�O�����g�̐���
  z_rate1 <- matrix(0, nrow=hhpt, ncol=seg_u)
  Zi1 <- matrix(0, nrow=hhpt, ncol=seg_u)
  z1_vec <- rep(0, hhpt)
  rf02 <- matrix(0, nrow=seg_u, ncol=seg_u)
  
  for(j in 1:max_pt){
    if(j==1){
      #�Z�O�����g�����m���𐄒�
      n <- n_time[j]
      LLs <- matrix(gamma1, nrow=n, ncol=seg_u, byrow=T) * LLi1[index_t11, ]   #�d�ݕt���ޓx
      matrix(gamma1, nrow=n, ncol=seg_u, byrow=T)
      z_rate1[index_t11, ] <- LLs / rowSums(LLs)   #�����m��
      
      #�������z����Z�O�����g�𐶐�
      Zi1[index_t11, ] <- rmnom(n, 1, z_rate1[index_t11, ])
      z1_vec[index_t11] <- as.numeric(Zi1[index_t11, ] %*% 1:seg_u)
      
      #�������̃p�����[�^���X�V
      rf01 <- colSums(Zi1[index_t11, ])
      
    } else {
      
      #�Z�O�����g�̊����m��
      index <- index_t22[[j]]
      n <- n_time[j]
      LLs <- gamma2[z1_vec[index_t21[[j]]], , drop=FALSE] * LLi1[index, , drop=FALSE]   #�d�ݕt���ޓx
      z_rate1[index, ] <- LLs / rowSums(LLs)   #�����m��
      
      #�������z���Z�O�����g�𐶐�
      Zi1[index, ] <- rmnom(n, 1, z_rate1[index, ])
      z1_vec[index] <- as.numeric(Zi1[index, ] %*% 1:seg_u)
      
      #�������̃p�����[�^���X�V
      rf02 <- rf02 + t(Zi1[index_t21[[j]], , drop=FALSE]) %*% Zi1[index, , drop=FALSE]   #�}���R�t����
    }
  }
  user_vec <- Zi1
  
  #�}���R�t���ڍs��̃p�����[�^���X�V
  rf1 <- rf01 + tau1
  rf2 <- rf02 + tau2
  gamma1 <- as.numeric(extraDistr::rdirichlet(1, rf1))
  gamma2 <- extraDistr::rdirichlet(seg_u, rf2)
  
  
  ##�A�C�e�����ƂɃZ�O�����g�����𐶐�
  #�Z�O�����g���Ƃɑ������z�̖ޓx���v�Z
  y2 <- sparse_data_T %*% user_vec   #���[�U�Z�O�����g���Ƃ̍w���p�x
  LLi2 <- as.matrix(y2 %*% log(t(oldphi)))
  
  #logsumexp�̖ޓx���v�Z
  LLi_max <- rowMaxs(LLi2)
  r_matrix <- matrix(omega1, nrow=item, ncol=seg_i, byrow=T)
  
  #�Z�O�����g�����m���̐���ƃZ�O�����g�̐���
  expl <- r_matrix * exp(LLi2 - LLi_max)
  z2_rate <- expl / rowSums(expl)   #�Z�O�����g�����m��
  item_vec <- Zi2 <- rmnom(item, 1, z2_rate)   #�������z����Z�O�����g�����𐶐�
  z2_vec <- as.numeric(Zi2 %*% 1:seg_i)

  #�������̍X�V
  rf3 <- colSums(Zi2) + tau3
  omega1 <- as.numeric(extraDistr::rdirichlet(1, rf3))
  
  
  ##�Z�O�����g�C���f�b�N�X���쐬
  index1 <- list()
  index2 <- list()
  for(i in 1:seg_u){index1[[i]] <- which(z1_vec==i)}
  for(j in 1:seg_i){index2[[j]] <- which(z2_vec==j)}

  
  ##���[�U�[����уA�C�e���̃p�����[�^���T���v�����O
  #�f�B�N�������z���烆�[�U�[�Z�O�����g�̃p�����[�^���T���v�����O
  freq_user <- matrix(0, nrow=seg_u, ncol=seg_i)
  for(i in 1:seg_u){
    x <- sparse_data[index1[[i]], , drop=FALSE]
    for(j in 1:seg_i){
      freq_user[i, j] <- sum(x[, index2[[j]]])
    }
  }
  freq_item <- t(freq_user)   #�A�C�e���̃p�����[�^�̓��[�U�[�p�����[�^�𔽓]�����邾��

  
  #�f�B�N�������z����A�C�e���Z�O�����g�̃p�����[�^���T���v�����O
  oldtheta <- extraDistr::rdirichlet(seg_u, freq_user + alpha1)
  oldphi <- extraDistr::rdirichlet(seg_i, freq_item + alpha2)
  
  ##�T���v�����O���ʂ̊i�[�ƕ\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[, , mkeep] <- oldtheta 
    PHI[, , mkeep] <- oldphi
    GAMMA1[mkeep, ] <- gamma1
    GAMMA2[, , mkeep] <- gamma2
    OMEGA1[mkeep, ] <- omega1
    
    if(rp >= 500){
    SEG1 <- SEG1 + Zi1 
    SEG2 <- SEG2 + Zi2 
    }
    
    if(rp%%disp==0){
      print(rp)
      print(round(rbind(gamma1, gamma01), 3))
      print(round(cbind(gamma2, gamma02), 3))
      print(round(cbind(oldtheta, thetat), 3))
      print(round(cbind(oldphi, phit), 3))
    }
  }
}
