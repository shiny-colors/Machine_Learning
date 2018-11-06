#####���Ӊ������f�B�N�����������f��#####
library(MASS)
library(vcd)
library(gtools)
library(caret)
library(reshape2)
library(plyr)
library(ggplot2)
library(lattice)

####�f�[�^�̔���####
##���f���̐ݒ�
k <- 5   #�Z�O�����g��
hh <- 500   #�Z�O�����g�̃T���v����
n <- hh * k   #���T���v����
w <- 200   #�Z�O�����g�ʂ̃p�����[�^��

#�Z�O�����g�̐ݒ�
seg.z <- rep(1:k, rep(hh, k))

##�l���Ƃ̕p�x��ݒ�
freq <- round(exp(rnorm(n, 4.75, 0.18)), 0)
summary(freq)
hist(freq, col="grey", main="���[�U�[���Ƃ̉{����")

##�Z�O�����g���ƂɊm���x�N�g�����`
#�f�B�N�������z���m���W����ݒ�
x <- rdirichlet(k, rep(0.1, w))
round(x, 3)   #�����������m�����m�F

##�����������m���ɂ��ƂÂ��p�x�f�[�^�𐶐�
#�������z���{�������𐶐�
Y <- matrix(0, nrow=n, ncol=w)
for(i in 1:n){
  Y[i, ] <- t(rmultinom(1, freq[i], x[seg.z[i], ]))
}
matrix(as.numeric(unlist(by(Y, seg.z, colSums))), nrow=k, ncol=w, byrow=T)   #�Z�O�����g���Ƃ̏o�������v�Z


####���Ӊ��M�u�X�T���v�����O�ő����f�B�N�������z�̃Z�O�����g�𐶐�####
##MCMC�A���S���Y���̐ݒ�
R <- 5000
keep <- 2
sbeta <- 1.5

##�n�C�p�[�p�����[�^�̐ݒ�
alpha <- 5
beta <- 5
beta.array <- array(beta, dim=c(k, w, n))


####���Ӊ��M�u�X�T���v�����O�ŃZ�O�����g�𐶐�####
##z�̏����l��ݒ�
#�����l�ˑ�����̂ŁA�ǂ������l��������܂Ŕ���������
for(t in 1:1000){
  print(t)
  #���[�U�[���Ƃ̃Z�O�����g�������s��`���ɕύX
  Z <- t(rmultinom(n, 1, rep(1/k, k)))
  z <- Z %*% 1:k
  
  zu <- matrix(0, nrow=n, ncol=k)
  za <- matrix(0, nrow=n, ncol=k)
  zw <- array(0, dim=c(k, w, n))
  
  for(i in 1:n){
    zu[i, z[i]] <- 1 
    za[i, z[i]] <- freq[i]
    zw[z[i], , i] <- Y[i, ]
  }
  
  ##�f�[�^�ƃA���S���Y������̊i�[�p�p�����[�^��ݒ�
  Nd <- matrix(freq, nrow=n, ncol=k)   #���[�U�[���Ƃ̕p�x
  
  #���[�U�[����уA�C�e�����Ƃ̕p�x
  Ndv <- array(0, dim=c(k, w, n))   
  for(i in 1:n){
    Ndv[, , i] <- matrix(Y[i, ], nrow=k, ncol=w, byrow=T)
  }
  Fdv <- ifelse(Ndv > 0, 1, 0)
  
  #���[�U�[���ƂɃA�C�e���̏o�����������擾
  index.n <- list()
  for(i in 1:n){
    index.n[[i]] <- subset(1:length(Y[i, ]), Y[i, ] > 0)
  }
  
  #���[�U�[���Ƃ̕p�x�s��
  freqM <- matrix(freq, nrow=n, ncol=k)
  
  #�p�x�̑������z��
  Y_array <- array(0, dim=c(k, w, n))
  for(i in 1:n){
    Y_array[, , i] <- matrix(Y[i, ], nrow=k, ncol=w, byrow=T)
  }
  
  #�|���A���z�̃p�������[�^�p
  lgamma3 <- matrix(0, nrow=n, ncol=k)
  lgamma4 <- matrix(0, nrow=n, ncol=k)
  
  ##�T���v�����O���ʂ̊i�[�p
  Seg_Z <- matrix(0, nrow=R/keep, ncol=n)
  Prob <- array(0, dim=c(n, k, R/keep))
  
  ##���v�ʂ̏����l���v�Z
  #���[�U�[�̃Z�O�����g����������
  Zu <- matrix(as.numeric(table(z)), nrow=n, ncol=k, byrow=T)
  Zul <- Zu - zu
  
  #���p�x�̃Z�O�����g����������
  Za <- matrix(as.numeric(by(Y, z, sum)), nrow=n, ncol=k, byrow=T)
  Zal <- Za - za
  
  #w�̃Z�O�����g����������
  Zw <- array(matrix(unlist(by(Y, z, colSums)), nrow=k, ncol=w, byrow=T), dim=c(k, w, n))
  Zwl <- Zw - zw
  
  
  ##��������mcmc�T���v�����O�����s
  for(rp in 1:R){
  
    ##�M�u�X�T���v�����O�̊m���̍X�V�����v�Z
    #��1���q���v�Z
    D_alpha <- log(Zul + alpha)
    
    #��2���q���v�Z
    lgamma1 <- lgamma(Zal + beta*w)
    lgamma2 <- lgamma(Zal + Nd + beta*w)
    
    #��3���q���v�Z
    lg4 <- Zwl + beta.array*Fdv
    lg3 <- lg4 + Ndv 
    
    for(i in 1:n){
      lgamma3[i, ] <- rowSums(lgamma(lg3[, index.n[[i]], i]))
      lgamma4[i, ] <- rowSums(lgamma(lg4[, index.n[[i]], i]))
    }
    
    ##�m���̌v�Z�ƃZ�O�����g�̐���
    #�Z�O�����g�����m�����v�Z
    #�ޓx�����������Ȃ��悤�ɒ萔��������
    lgamma_c <- D_alpha + lgamma1 - lgamma2 + lgamma3 - lgamma4
    lgamma_mean <- matrix(apply(lgamma_c, 1, mean), nrow=n, ncol=k)
    lgamma <- exp(lgamma_c - lgamma_mean)
    lgamma[is.infinite(lgamma) & lgamma > 0] <- 10^300
    lgamma[is.infinite(lgamma) & lgamma < 0] <- -10^300
    
    #�m���̌v�Z
    Pr <- lgamma / rowSums(lgamma)
  
    #�������z���Z�O�����g�𐶐�
    Z <- t(apply(Pr, 1, function(x) rmultinom(1, 1, x)))
    z <- Z %*% 1:k
    
    if(length(table(z)) < k) {break}   #�Z�O�����g���k�ނ�����break����
    
    ##���v�ʂ��X�V
    #�Z�O�����g����������
    #���[�U�[�̃Z�O�����g����������
    Zu <- matrix(as.numeric(table(z)), nrow=n, ncol=k, byrow=T)
    zu <- Z
    Zul <- Zu - zu
    
    #���p�x�̃Z�O�����g����������
    Za <- matrix(as.numeric(by(Y, z, sum)), nrow=n, ncol=k, byrow=T)
    za <- Z * freqM   
    Zal <- Za - za
    
    #w�̃Z�O�����g����������
    Zw <- array(matrix(unlist(by(Y, z, colSums)), nrow=k, ncol=w, byrow=T), dim=c(k, w, n))
    Z_array <- array(as.numeric(t(Z[rep(1:n, rep(w, n)), ])), dim=c(k, w, n))   
    zw <- Z_array * Y_array
    Zwl <- Zw - zw
    
    ##�p�����[�^���i�[
    if(rp%%keep==0){
      print(rp)
      mkeep <- rp/keep
      Seg_Z[mkeep, ] <- z
      Prob[, , mkeep] <- Pr
      print(table(z))
    }
  }
  if(rp==R) {break}   #�ݒ肵���T���v�����O���ɒB���Ă���΁Amcmc���[�v���I������
}

####���茋�ʂ̗v��ƓK���x�̊m�F####
burnin <- 2000/keep   #�o�[���C������

#�Z�O�����g�̃T���v�����O����
Seg_Z[burnin:(R/keep), c(1:4, 501:504, 1001:1004, 1501:1504, 2001:2004)]

##�m���̃T���v�����O����
round(colMeans(t(Prob[1, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[501, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[1001, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[1501, , burnin:(R/keep)])), 3)
round(colMeans(t(Prob[2001, , burnin:(R/keep)])), 3)