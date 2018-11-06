#####��������-�񍀕��z���f��#####
library(MASS)
library(gtools)
library(reshape2)
library(plyr)
library(dplyr)
library(knitr)
library(ggplot2)


####�f�[�^�̔���####
##���f���̐ݒ�
k <- 5   #�Z�O�����g��
n <- 500   #�Z�O�����g�̃T���v����
N <- n * k   #���T���v����
seg <- rep(1:k, rep(n, k))   #�Z�O�����g����
ns1 <- rpois(N, 30)   #�������z�̌l���Ƃ̕p�x
ns2 <- rpois(N, 40)   #�񍀕��z�̌l���Ƃ̕p�x
th <- 10   #�Z�O�����g�ʂ̃p�����[�^��


##�m���x�N�g�����`���āA�����ϐ��𔭐�
P1 <- matrix(0, nrow=k, ncol=th) 
P2 <- rep(0, k)
X1.list <- list()
X2 <- c()

for(j in 1:k){
  r <- ((j-1)*n+1):((j-1)*n+n)
  
  #�m�����v�Z
  p <- rgamma(th, 0.8, 0.2)
  P1[j, ] <- p / sum(p)
  P2[j] <- runif(1, 0.15, 0.7)
  
  #�����ϐ��𔭐�
  X1.list[[j]] <- t(apply(cbind(ns1[r], 1), 1, function(x) rmultinom(1, x[1], P1[j, ])))
  X2 <- c(X2, apply(cbind(ns2[r], 1), 1, function(x) rbinom(1, x[1], P2[j])))
}

X1 <- do.call(rbind, X1.list)


##�f�[�^���������Č��ʂ��W�v
#�����f�[�^�̗v��
by(X1, seg, function(x) round(colMeans(x), 3))
by(X1, seg, function(x) summary(x))
by(X1, seg, function(x) round(colSums(x)/sum(x), 3))

#�񍀃f�[�^�̗v��
by(X2, seg, function(x) round(mean(x), 3))
by(X2, seg, function(x) summary(x))
by(cbind(X2, ns2), seg, function(x) round(sum(x[1])/sum(x[2]), 3))


####EM�A���S���Y���ō�������-�񍀕��z���f���𐄒�####
##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta1, theta2, n1, n2, X1, X2, r, k){
  LLind <- matrix(0, nrow=nrow(X1), ncol=k)
  for(i in 1:k){
    Li1 <- apply(cbind(n1, X1), 1, function(x) dmultinom(x[-1], x[1], theta1[i, ]))   #�������z�̖ޓx
    Li2 <- apply(cbind(n2, X2), 1, function(x) dbinom(x[-1], x[1], theta2[i]))   #�񍀕��z�̖ޓx
    Li <- Li1 * Li2
    LLind[, i] <- Li
  }
  
  LLho <- matrix(r, nrow=nrow(X1), ncol=k, byrow=T) * LLind   #�ϑ��f�[�^�̖ޓx
  z <- LLho / matrix(apply(LLho, 1, sum), nrow=nrow(X1), ncol=k)   #���ݕϐ�z�̌v�Z
  LLosum <- sum(log(apply(matrix(r, nrow=nrow(X1), ncol=k, byrow=T) * LLind, 1, sum)))   #�ϑ��f�[�^�̑ΐ��ޓx�̌v�Z
  rval <- list(LLob=LLosum, z=z, LL=LLind, Li1=Li1, Li2=Li2)
  return(rval)
}

#�����l�̐ݒ�
iter <- 0
k <- 5   #�Z�O�����g��

##theta�̏����l�̐ݒ�
theta1 <- matrix(0, nrow=k, ncol=th)
theta2 <- c()

for(j in 1:k){
  p <- runif(th, 0.1, 1.5)
  theta1[j, ] <- p / sum(p)
  theta2 <- c(theta2, runif(1, 0.25, 0.7))
}

##������r�̏����l
r <- c(0.15, 0.15, 0.25, 0.2, 0.25)

#�ΐ��ޓx�̏�����
L <- LLobz(theta1=theta1, theta2=theta2, n1=ns1, n2=ns2, X1=X1, X2=X2, r=r, k=k)
LL1 <- L$LLob
z <- L$z
round(z, 3)

#�X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 1  

##EM�A���S���Y��
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�
  #E�X�e�b�v�̌v�Z
  z <- L$z   #���ݕϐ�z�̏o��
  
  #M�X�e�b�v�̌v�Z�ƍœK��
  #theta�̐���
  theta <- matrix(0, nrow=k, ncol=th)
  for(j in 1:k){
    #���S�f�[�^�̑ΐ��ޓx����theta�̐���ʂ��v�Z
    thetaseg1 <- colSums(matrix(z[, j], nrow=nrow(X1), ncol=th) * X1) / sum(z[, j] * ns1)   #�d�ݕt���������z�̍Ŗސ���
    thetaseg2 <- sum(z[, j]*X2) / sum(z[, j]*ns2)   #�d�ݕt���񍀕��z�̍Ŗސ���
    
    theta1[j, ] <-thetaseg1
    theta2[j] <- thetaseg2  
  }
  
  #�������𐄒�
  r <- apply(z, 2, sum) / nrow(X1)
  
  #�ϑ��f�[�^�̑ΐ��ޓx���v�Z
  L <- LLobz(theta1=theta1, theta2=theta2, n1=ns1, n2=ns2, X1=X1, X2=X2, r=r, k=k)
  LL <- L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
  iter <- iter+1   
  dl <- LL-LL1
  LL1 <- LL
  print(LL)
}


####��������-�񍀕��z���f���̐��茋��####
round(theta1, 3)   #theta1�̐����
round(P1, 3)   #theta1�̐^�̒l
round(theta2, 3)   #theta1�̐����
round(P2, 3)   #theta1�̐^�̒l
round(r, 3)   #�������̐����
round(z, 3)   #�l�ʂ̃Z�O�����g�ւ̏����m��

L$LLob   #�ϑ��f�[�^�̑ΐ��ޓx
-2*(L$LLob) + 2*(k*nrow(theta1)+length(theta2) + 1)   #AIC



