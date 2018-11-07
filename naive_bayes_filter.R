#####�i�C�[�u�x�C�Y�t�B���^#####
library(MASS)
library(e1071)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
Nnm <- 1500   #����f�[�^�̃T���v����
Nab <- 1500   #�ُ�f�[�^�̃T���v����
N <- 3000   #���T���v����
k <- 200   #��b��
w <- 500   #����������̒P�ꐔ
z <- c(rep(0, Nnm), rep(1, Nab))   #����ƈُ�̃C���f�B�P�[�^�[

#���ނɊ֘A���Ȃ��f�[�^�̕��z�̔���
W1 <- matrix(0, N, 3/4*k) 
for(i in 1:(3/4*k)){
  p1 <- runif(1, 0.01, 0.25)
  W1[, i] <- rbinom(N, 1, p1)
}

#����f�[�^�Ɋ֘A����f�[�^�̔���
Wnm <- matrix(0, N, 1/4*k/2)
for(i in 1:(1/4*k/2)){
  pnm1 <- runif(1, 0.075, 0.3)
  pnm2 <- runif(1, 0.005, 0.1)
  Wnm1 <- rbinom((N/2), 1, pnm1)
  Wnm2 <- rbinom((N/2), 1, pnm2)
  Wnm.r <- c(Wnm1, Wnm2)
  Wnm[, i] <- Wnm.r
}

#�ُ�f�[�^�Ɋ֘A����f�[�^�̔���
Wab <- matrix(0, N, 1/4*k/2)
for(i in 1:(1/4*k/2)){
  pab1 <- runif(1, 0.01, 0.1)
  pab2 <- runif(1, 0.15, 0.3)
  Wab1 <- rbinom((N/2), 1, pab1)
  Wab2 <- rbinom((N/2), 1, pab2)
  Wab.r <- c(Wab1, Wab2)
  Wab[, i] <- Wab.r
}

#�f�[�^�̌���
Wz <- cbind(z, W1, Wnm, Wab)

#�f�[�^�̏W�v
by(Wnm, z, colSums)
by(Wab, z, colSums)


####�i�C�[�u�x�C�Y���ފ�####
#�w�K�f�[�^�ƃe�X�g�f�[�^�ɕ�����
Wz.learn <- rbind(Wz[1:1000, ], Wz[1501:2500, ])
Wz.test <- rbind(Wz[1001:1500, ], Wz[2501:3000, ])
z.learn <- c(z[1:1000], z[1501:2500])
z.test <- c(z[1001:1500], z[2501:3000])

##�i�C�[�u�x�C�Y�̖ޓx�֐�
#�P��p�x�𐔂���
wnm <- colSums(Wz.learn[z.learn==0, -1])   #����f�[�^�ł̒P��p�x
wab <- colSums(Wz.learn[z.learn==1, -1])   #�ُ�f�[�^�ł̒P��p�x

#�����t���m��
round(pnm <- wnm / (wnm + wab), 3)   #����f�[�^�̏����t���m��
round(pab <- wab / (wnm + wab), 3)   #�ُ�f�[�^�̏����t���m��

#���O�m���@
prior <- 0.5

##�x�C�Y�t�B���^�Ő���f�[�^�A�ُ�f�[�^�𔻒�(���O�m���͓����@)
score <- c()
for(i in 1:nrow(Wz.test)){
  #����A�ُ�̖ޓx
  post.nm <- Wz.test[i, -1]*pnm
  post.ab <- Wz.test[i, -1]*pab
  
  #�ޓx���[���ȊO�̃f�[�^�����o��
  post.nmz <- subset(post.nm, post.nm > 0)
  post.abz <- subset(post.ab, post.ab > 0)
  
  #�X�R�A���v�Z
  score <- c(score, sum(log(post.abz))-sum(log(post.nmz)))
}

#����
round(score, 3)
(res <- ifelse(score > 0, 1, 0))
(res.table <- table(z.test, res))   #�딻�ʕ\
c(res.table[1, 1] / sum(res.table[1, ]), res.table[2, 2] / sum(res.table[2, ]))   #���ޕʂ̐�����
sum(diag(res.table)) / sum(res.table)   #�S�̂ł̐�����

