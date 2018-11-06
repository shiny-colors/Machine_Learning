#####�����������z���f��(���O�T��exp���p)#####
library(MASS)
library(vcd)
library(gtools)
library(matrixStats)
library(Matrix)
library(extraDistr)
library(caret)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(lattice)

#set.seed(6439)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 20000   #���[�U�[��
k <- 10   #�Z�O�����g��
v <- 400   #�A�C�e��
s <- rpois(hh, 20)   #�w����
s[s < 5] <- ceiling(runif(sum(s < 5), 5, 10))
f <- sum(s)

##�p�����[�^�̐ݒ�
theta <- thetat <- extraDistr::rdirichlet(1, rep(10, k))
phi <- phit <- extraDistr::rdirichlet(k, rep(0.3, v))

##�������z���f�[�^�𐶐�
Z_list <- list()
Data <- matrix(0, nrow=hh, ncol=v)

for(i in 1:hh){
  #�Z�O�����g�����𐶐�
  z <- rmnom(1, 1, theta)
  z_vec <- as.numeric(z %*% 1:k)
  Z_list[[i]] <- z
  
  #�w���f�[�^�𐶐�
  w <- rmnom(s[i], 1, phi[z_vec, ])
  Data[i, ] <- colSums(w)
}
Z <- do.call(rbind, Z_list)

####EM�A���S���Y���ō����������z�𐄒�####
const <- lfactorial(s) - rowSums(lfactorial(Data))   #�������z�̖��x�֐��̑ΐ��ޓx�̒萔
r <- colMeans(Z)

##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(Data, phi, r, const, hh, k){
  
  #�������z�̑ΐ��ޓx
  log_phi <- log(t(phi))
  LLi <- const + Data %*% log_phi
  
  #logsumexp�̖ޓx
  LLi_max <- matrix(apply(LLi, 1, max), nrow=hh, ncol=k)
  r_matrix <- matrix(r, nrow=hh, ncol=k, byrow=T)
  
  #�����m���̃p�����[�^��ݒ�
  expl <- r_matrix * exp(LLi - LLi_max)
  expl_log <- log(expl)
  expl_max <- matrix(log(max(expl[1, ])), nrow=hh, ncol=k)
  z <- exp(expl_log - (log(rowSums(exp(expl_log - expl_max))) + expl_max))   #�Z�O�����g�����m��
  
  #�ϑ��f�[�^�̑ΐ��ޓx
  r_log <- matrix(log(r), nrow=hh, ncol=k, byrow=T)
  LLosum <- sum(log(rowSums(exp(r_log + LLi))))   #�ϑ��f�[�^�̑ΐ��ޓx
  rval <- list(LLob=LLosum, z=z, LL=LLi)
  return(rval)
}

##�p�����[�^�̏����l]
#phi�̏����l
alpha0 <- colSums(Data) / sum(Data)
phi <- extraDistr::rdirichlet(k, alpha0*1000)

#�������̏����l
r <- rep(1/k, k)

#�ϑ��f�[�^�̑ΐ��ޓx�̏�����
L <- LLobz(Data, phi, r, const, hh, k)
LL1 <- L$LLob
z <- L$z

#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.5
iter <- 0 

##EM�A���S���Y���őΐ��ޓx���ő剻
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  
  #M�X�e�b�v�̌v�Z�ƍœK��
  #phi�̐���
  phi <- matrix(0, nrow=k, ncol=v)
  for(j in 1:k){
    #���S�f�[�^�̑ΐ��ޓx����phi�̐���ʂ��v�Z
    phi[j, ] <- colSums(matrix(z[, j], nrow=hh, ncol=v) * Data) / sum(z[, j] * s)   #�d�ݕt���������z�̍Ŗސ���
  }
  
  #�������𐄒�
  r <- apply(z, 2, sum) / hh
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  phi[phi==0] <- min(phi[phi > 0])
  L <- LLobz(Data, phi, r, const, hh, k)
  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}

####�����������z���f���̐��茋��####
##EM�A���S���Y���̐��茋�ʂƐ^�l�̔�r
round(cbind(t(phi), t(phit)), 3)   #phi�̐���l��phi�̐^�l
round(cbind(Z %*% 1:k, apply(z, 1, which.max), z), 3)   #���ݕϐ�z�ƃZ�O�����g����
round(r, 3)   #������

##�K���x�̊m�F
L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
-2*(L$LLob) + 2*length(phi)   #AIC


