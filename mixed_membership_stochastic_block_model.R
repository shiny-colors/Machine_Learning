#####���������o�V�b�v�m���I�u���b�N���f��#####
library(MASS)
library(bayesm)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(reshape2)
library(qrmtools)
library(slfm)
library(caret)
library(dplyr)
library(foreach)
library(ggplot2)
library(lattice)

set.seed(4965723)

####�f�[�^�̔���####
k <- 8   #������
d <- 2000   #�m�[�h��

#�m�[�h�\����ݒ�
theta <- extraDistr::rdirichlet(d, rep(0.2, 2*k))
omega <- rbeta(d, 2.5, 5.0)
gamma <- rbeta(2*k, 1.2, 1.5)
index_list <- list()
id_list <- list()

for(i in 1:d){
  #�m�[�h�m����ݒ�
  par <- matrix(theta[i, ], nrow=d-1, ncol=2*k) * theta[-i, ]
  z_rate <- par / rowSums(par)
  z_vec <- as.numeric(rmnom(d-1, 1, z_rate) %*% 1:(2*k))
  
  #�m�[�h�𐶐�
  index <- (1:d)[-i]
  flag <- rbinom(d-1, 1, omega[i] * gamma[z_vec])
  id_list[[i]] <- rep(i, d-1)
  index_list[[i]] <- index * flag
}

#���X�g��ϊ�
node_vec0 <- unlist(index_list)
index_node <- which(node_vec0 > 0)
node_vec <- node_vec0[index_node]
id_vec <- unlist(id_list)[index_node]

#�C���f�b�N�X���쐬
id_list1 <- list()
d_vec1 <- list()
id_list2 <- list()
d_vec2 <- list()
w <- rep(0, d)
for(i in 1:d){
  id_list1[[i]] <- which(id_vec==i)
  id_list2[[i]] <- which(node_vec==i)
  d_vec1[[i]] <- rep(1, length(id_list1[[i]]))
  d_vec2[[i]] <- rep(1, length(id_list2[[i]]))
  w[i] <- length(id_list1[[i]])
}
f <- sum(w)


##���f���Ɋ�Â��f�[�^�𐶐�
#�p�����[�^�̎��O���z
k0 <- k/2
alpha1 <- rep(0.15, k)
beta01 <- 5.0
beta02 <- 5.5
beta03 <- 1.5
beta04 <- 5.0


#�p�����[�^�𐶐�
phi <- matrix(0, nrow=k, ncol=k)
for(j in 1:k){
  index_phi <- cbind(rep(j, k), 1:k)
  for(l in 1:k){
    #�x�[�^���z���烊���N�m���𐶐�
    index <- index_phi[l, ]
    if(index[1]==index[2]){
      phi[index[1], index[2]] <- rbeta(1, beta01, beta02)
    } else {
      phi1 <- rbeta(1, beta03, beta04)     
      phi2 <- rbeta(1, (phi1)*k0, (1-phi1)*k0)
      phi[index[1], index[2]] <- phi1; phi[index[2], index[1]] <- phi2
    }
  }
}
phit <- phi
theta1 <- thetat1 <- extraDistr::rdirichlet(d, alpha1)
theta2 <- thetat2 <- extraDistr::rdirichlet(d, theta1 + 0.05)

#�m�[�h���ƂɃ����N�𐶐�
y <- rep(0, f)
Z1_list <- list()
Z2_list <- list()

for(i in 1:d){
  #�g�s�b�N�𐶐�
  z1 <- rmnom(w[i], 1, theta1[i, ])
  z2 <- rmnom(w[i], 1, theta2[node_vec[id_list1[[i]]], ])
  z1_vec <- as.numeric(z1 %*% 1:k)
  z2_vec <- as.numeric(z2 %*% 1:k)
  
  #�����N�𐶐�
  phi_vec <- rowSums(phi[z1_vec, ] * z2)
  y[id_list1[[i]]] <- rbinom(w[i], 1, phi_vec)
  
  #�f�[�^���i�[
  Z1_list[[i]] <- z1
  Z2_list[[i]] <- z2
}

#���X�g��ϊ�
Z1 <- do.call(rbind, Z1_list)
Z2 <- do.call(rbind, Z2_list)
index_link <- which(y==1)


####�}���R�t�A�������e�J�����@��MMSB�𐄒�####
##�A���S���Y���̐ݒ�
R <- 5000
keep <- 2  
iter <- 0
burnin <- 1000/keep
disp <- 10

##���O���z�̐ݒ�
alpha1 <- k
beta1 <- 0.1
beta2 <- 0.1

##�p�����[�^�̐^�l
theta1 <- thetat1
theta2 <- thetat2
phi <- phit

##�p�����[�^�̏����l
theta1 <- extraDistr::rdirichlet(d, rep(5.0, k))
theta2 <- extraDistr::rdirichlet(d, rep(5.0, k))
phi <- matrix(rbeta(k*k, 1.0, 2.5), nrow=k, ncol=k)
diag(phi) <- 0.5
Zi1 <- rmnom(f, 1, theta1[id_vec, ])
Zi2 <- rmnom(f, 1, theta2[node_vec, ]) 
z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)

##�p�����[�^�̊i�[�p�z��
THETA1 <- array(0, dim=c(d, k, R/keep))
THETA2 <- array(0, dim=c(d, k, R/keep))
PHI <- array(0, dim=c(d, k, R/keep))
SEG1 <- matrix(0, nrow=f, ncol=k)
SEG2 <- matrix(0, nrow=f, ncol=k)


##�C���f�b�N�X���쐬
index10 <- id_vec[-index_link]
index11 <- id_vec[index_link]
index20 <- node_vec[-index_link]
index21 <- node_vec[index_link]

##�ΐ��ޓx�̊�l
phi_mu <- mean(y)
LLst <- sum(dbinom(y, 1, phi_mu, log=TRUE))


####�M�u�X�T���v�����O�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�m�[�h�̃g�s�b�N�𐶐�
  #�g�s�b�N���Ƃ̊��Җޓx���v�Z
  Li1 <- matrix(0, nrow=f, ncol=k)
  Li2 <- matrix(0, nrow=f, ncol=k)
  
  for(j in 1:k){
    Li1[-index_link, j] <- theta1[index10, j] * (1-phi)[j, ][z2_vec[index20]]
    Li1[index_link, j] <- theta1[index11, j] * phi[j, ][z2_vec[index21]]
    Li2[-index_link, j] <- theta2[index20, j] * (1-phi)[, j][z1_vec[index10]]
    Li2[index_link, j] <- theta2[index21, j] * phi[, j][z1_vec[index11]]
  }

  #�������z����m�[�h�̃g�s�b�N�𐶐�
  z1_rate <- Li1 / rowSums(Li1); z2_rate <-  Li2 / rowSums(Li2)   #���ݕϐ�z
  Zi1 <- rmnom(f, 1, z1_rate); Zi2 <- rmnom(f, 1, z2_rate)   #�g�s�b�N���T���v�����O
  z1_vec <- as.numeric(Zi1 %*% 1:k); z2_vec <- as.numeric(Zi2 %*% 1:k)
  Zi1_T <- t(Zi1); Zi2_T <- t(Zi2)
  
  ##�p�����[�^���T���v�����O
  #�g�s�b�N���z���T���v�����O
  wsum01 <- wsum02 <- matrix(0, nrow=d, ncol=k)
  for(i in 1:d){
    wsum01[i, ] <- Zi1_T[, id_list1[[i]]] %*% d_vec1[[i]] 
    wsum02[i, ] <- Zi2_T[, id_list2[[i]]] %*% d_vec2[[i]]
  }
  wsum1 <- wsum01 + alpha1; wsum2 <- wsum02 + alpha1
  theta1 <- extraDistr::rdirichlet(d, wsum1)   #�f�B���N�����z����theta���T���v�����O
  theta2 <- extraDistr::rdirichlet(d, wsum2)
  
  #�����N�m�����T���v�����O
  vsum0 <- t((1-y) * Zi1) %*% ((1-y) * Zi2) + beta1
  vsum1 <- t(y * Zi1) %*% (y * Zi2) + beta2
  phi <- matrix(rbeta(k*k, vsum1, vsum0), nrow=k, ncol=k)
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    THETA1[, , mkeep] <- theta1
    THETA2[, , mkeep] <- theta2
    PHI[, , mkeep] <- phi
  }
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG2 <- SEG2 + Zi2
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    phi_vec <- rep(0, f)
    z1_vec <- as.numeric(Zi1 %*% 1:k)
    z2_vec <- as.numeric(Zi2 %*% 1:k)
    for(i in 1:f){
      phi_vec[i] <- phi[z1_vec[i], z2_vec[i]]
    }
    #�T���v�����O���ʂ̕\��
    print(rp)
    print(c(sum(dbinom(y, 1, phi_vec, log=TRUE)), LLst))
    print(round(cbind(phi, phit), 3))
  }
}


####�T���v�����O���ʂ̉����Ɨv��####
#�g�s�b�N���z�̃T���v�����O����
matplot(t(THETA1[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA1[5, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[10, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(THETA2[15, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

#�����N�m���̃p�����[�^
matplot(t(PHI[1, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[2, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[3, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(PHI[4, , ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
