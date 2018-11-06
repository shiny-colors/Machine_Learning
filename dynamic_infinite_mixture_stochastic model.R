#####�x�C�W�A���������ݐ��ړ񍀕��z���f��#####
library(MASS)
library(lda)
library(RMeCab)
library(bayesm)
library(MCMCpack)
library(extraDistr)
library(matrixStats)
library(gtools)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(45327)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 5000   #���[�U�[��
pt <- 7   #�ϑ����Ԑ���
hhpt <- hh*pt   #���T���v����
seg <- 7   #�Z�O�����g��
k <- 50   #�ϑ��ϐ���

##ID�̐ݒ�
id <- rep(1:hh, rep(pt, hh))
time <- rep(1:pt, hh)
ID <- data.frame(no=1:hhpt, id, time)

##�Z�O�����g�̐ݒ�
#�Z�O�����g���ڊm���s��̐ݒ�
P_seg0 <- matrix(0, nrow=seg, ncol=seg)
for(i in 1:seg){
  rand <- runif(seg, 0.1, 4)
  P_seg0[i, ] <- rand
}
diag(P_seg0) <- runif(seg, 5, 25)   #�Ίp�s���u������
P_seg <- P_seg0 / matrix(rowSums(P_seg0), nrow=seg, ncol=seg)   #�m���ɒu������

#�Z�O�����g�𔭐�������
#����1�̃Z�O�����g��ݒ�
seg_m <- matrix(0, nrow=hh, ncol=pt)
seg_m[, 1] <- rep(1:seg, rep(ceiling(hh/seg), seg))[1:hh]

#����1~7�܂Œ����I�ɃZ�O�����g�𔭐�������
for(j in 2:pt){
  for(i in 1:hh){
    seg_m[i, j] <- extraDistr::rmnom(1, 1, P_seg[seg_m[i, j-1], ]) %*% 1:seg
  }
}
seg_vec <- as.numeric(t(seg_m))   #�Z�O�����g���x�N�g���ɕύX
r_rate0 <- do.call(rbind, tapply(seg_vec, ID$time, function(x) table(x)/sum(table(x))))   #�������̐^�l
r_rate0 <- matrix(as.numeric(r_rate0), nrow=pt, ncol=seg)

####�Z�O�����g�����Ɋ�Â������ϐ��𔭐�####
##�Z�O�����g�ʂ̃p�����[�^��ݒ�
Pr0 <- matrix(0, nrow=seg, ncol=k)
for(i in 1:seg){
  for(j in 1:k){
    tau1 <- runif(1, 0.5, 4.5)
    tau2 <- runif(1, 2.0, 9.5)
    Pr0[i, j] <- rbeta(1, tau1, tau2)
  }
}


##�񍀕��z���牞���ϐ��𔭐�
Data <- matrix(0, nrow=hhpt, ncol=k)
for(i in 1:seg){
  Data[seg_vec==i, ] <- matrix(rbinom(sum(seg_vec==i)*k, 1, Pr0[i, ]), nrow=sum(seg_vec==i), ncol=k, byrow=T)
}

#�f�[�^�̊m�F�ƏW�v
colSums(Data); colMeans(Data)
by(Data, seg_vec, function(x) round(colMeans(x), 3))
by(Data, seg_vec, function(x) round(colSums(x), 3))
by(Data, seg_vec, function(x) summary(x))


####�}���R�t�A�������e�J�����@�Ŗ������ݐ��ڊm�����f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 20000
keep <- 4 
sbeta <- 1.5
iter <- 0

##���O���z�̐ݒ�
tau <- c(1, 1)   #�x�[�^���z�̎��O���z
alpha <- 1   #CRP�̎��O���z

##�����l�̐ݒ�
seg0 <- 2   #�����Z�O�����g��2��
r <- c(0.5, 0.5)   #�������̏����l
par_mean <- colMeans(Data)   #CRP�p�̃p�����[�^

#�����Z�O�����g�̐ݒ�
z <- matrix(0, nrow=hhpt, ncol=seg0)
z0 <- rbinom(hhpt, 1, 0.5) + 1
for(i in 1:seg0){z[z0==i, i] <- 1}
z_vec <- z %*% 1:seg0

#�Z�O�����g���Ƃ̃p�����[�^
oldtheta <- matrix(rbeta(k*seg0, 2.0, 2.0), nrow=seg0, ncol=k)

##�p�����[�^�̊i�[�p
max_seg <- 20
Z <- matrix(0, nrow=hhpt, ncol=max_seg)
P <- array(0, dim=c(max_seg, k, R/keep))
storage.mode(Z) <- "integer"


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  LLind0 <- matrix(0, nrow=hhpt, ncol=nrow(oldtheta))
  for(j in 1:ncol(LLind0)){
    Li <- Data %*% log(oldtheta[j, ]) + (1-Data) %*% log(1-oldtheta[j, ])
    LLind0[, j] <- Li
  }
  
  #�V�������ݕϐ��̖ޓx�̌v�Z�Ɩޓx�̌���
  par <- colMeans(Data)
  LL_new <- Data %*% log(par) + (1-Data) %*% log(1-par)
  LLi0 <- cbind(LLind0, LL_new)
  LLi <- exp(LLi0 - max(LLi0))   #�ޓx�ɕϊ�
  
  ##���ݕϐ��̊����m���̎��O���z��CRP�̌v�Z
  r0 <- list()
  r <- matrix(0, nrow=hhpt, ncol=ncol(z)+1)
  
  #1���ł͑S�̂̍�������2���ȍ~�͑O���̍��������̗p
  for(j in 1:pt){
    if(j==1){
      index <- which(ID$time==j)
      r0[[j]] <- cbind(matrix(colSums(z), nrow=length(index), ncol=ncol(z), byrow=T) - z[index, ], alpha)
      r[index, ] <- r0[[j]] / (hhpt-1-alpha/pt)
    }
    if(j!=1){
      index_now <- which(ID$time==j)
      index_obs <- which(ID$time==j-1)
      r0[[j]] <- cbind(matrix(colSums(z[index_obs, ]), nrow=length(index_obs), ncol=ncol(z), byrow=T) - z[index_obs, ], alpha/pt)
      r[index_now, ] <- r0[[j]] / (length(index_obs)-1-alpha/pt)
    }
  }
  
  ##���ݕϐ��̊����m���̌v�Z�Ɛ��ݕϐ��̃T���v�����O
  gamma <- LLi * r
  z_rate <- gamma / rowSums(gamma)   #���ݕϐ�z�̊����m��
  z <- rmnom(hhpt, 1, z_rate)   #�������z������ݕϐ�z���T���v�����O
  z[is.nan(z)] <- 0
  z <- z[, colSums(z) > 2]
  
  ##�񍀕��z�̃p�����[�^���X�V
  oldtheta <- matrix(0, nrow=ncol(z), ncol=k) 
  for(j in 1:ncol(z)){
    
    #�p�����[�^�X�V�ɕK�v�ȃf�[�^�𒊏o
    index <- which(z[, j]==1)
    x <- Data[index, ]   #�����Z�O�����g�ɊY������f�[�^�𒊏o
    
    if(length(x) > k){
      y <- colSums(x)
      n <- nrow(x)
    } else {
      y <- x
      n <- 1
    }
    
    #�x�[�^���z����p�����[�^�𔭐�
    phi1 <- tau[1] + y
    phi2 <- tau[2] + n - y
    oldtheta[j, ] <- rbeta(k, phi1, phi2)
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    if(rp >= R/4){Z[, 1:ncol(z)] <- Z[, 1:ncol(z)] + z}   #�J��Ԃ������ő唽������1/4�𒴂�����p�����[�^���i�[
    P[1:nrow(oldtheta), , mkeep] <- oldtheta
    
    print(rp)
    print(round(colMeans(z), 3))
    print(round(rbind(oldtheta, Pr0)[, 1:15], 3))
  }
}


####�T���v�����O���ʂ̗v��Ɖ���####
burnin <- 2500/keep
RS <- R/keep

##�T���v�����O���ʂ̉���
matplot(t(P[1:seg, 1, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(P[1:seg, 2, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(P[1:seg, 3, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(P[1:seg, 4, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")
matplot(t(P[1:seg, 5, ]), type="l", xlab="�T���v�����O��", ylab="�p�����[�^")

##�T���v�����O���ʂ̗v��
round(rbind(apply(P[1:seg, , burnin:RS], c(1, 2), mean), Pr0), 3)   #���㕽��
round(apply(P[1:seg, , burnin:RS], c(1, 2), sd), 3)   #����W���΍�
round(cbind(Z[, colSums(Z) > 0] / rowSums(Z), seg=seg_vec), 3)   #���ݕϐ��̊����̎��㕪�z

