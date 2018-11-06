#####���ݐ��ڑ������z���f��#####
library(MASS)
library(mclust)
library(flexmix)
library(reshape2)
library(dplyr)
library(lattice)
library(ggplot2)
library(knitr)

#set.seed(573298)

####�f�[�^�̔���####
n <- 3000   #�T���v����
seg <- 5   #�Z�O�����g��
k <- 25   #�ϑ��ϐ���
period <- 5   #�ϑ�����
N <- n*period   #���T���v����
freq <- rpois(n*seg, 35)   #�ϑ���


##ID�̐ݒ�
id <- rep(1:n, rep(period, n))
time <- rep(1:period, n)
ID <- data.frame(no=1:(n*period), id=id, time=time)

##�Z�O�����g�̐ݒ�
#�Z�O�����g���ڊm���s��̐ݒ�
P_seg0 <- matrix(0, nrow=seg, ncol=seg)
for(i in 1:seg){
  rand <- runif(seg, 0.1, 3.5)
  P_seg0[i, ] <- rand
}
diag(P_seg0) <- runif(seg, 3.5, 25.0)   #�Ίp�s���u������
P_seg <- P_seg0 / matrix(rowSums(P_seg0), nrow=seg, ncol=seg)   #�m���ɒu������

#�Z�O�����g�𔭐�������
#����1�̃Z�O�����g�̐ݒ�
seg_m <- matrix(0, nrow=n, ncol=seg)
seg_m[, 1] <- rep(1:seg, rep(n/seg, seg))

#����1�`5�܂Œ����I�ɃZ�O�����g�𔭐�������
for(j in 2:seg){
  for(i in 1:n){
    seg_m[i, j] <- t(rmultinom(1, 1, P_seg[seg_m[i, j-1], ])) %*% 1:seg
  }
}
seg_vec <- as.numeric(t(seg_m))   #�Z�O�����g���x�N�g���ɕύX
table(seg_vec)
r_rate0 <- do.call(rbind, tapply(seg_vec, ID$time, function(x) table(x)/sum(table(x))))   #������
r_rate0 <- matrix(as.numeric(r_rate0), nrow=period, ncol=seg)


####�Z�O�����g�����Ɋ�Â������ϐ��𔭐�####
##�Z�O�����g�ʂ̊m���𔭐�
P <- matrix(0, nrow=seg, ncol=k)
for(i in 1:seg){
  p <- rgamma(k, 0.8, 0.2)
  P[i, ] <- p / sum(p)
}

##�������z���牞���ϐ��𔭐�
Y <- t(apply(cbind(freq, seg_vec), 1, function(x) rmultinom(1, x[1], P[x[2], ])))

#�f�[�^�̊m�F�ƏW�v
round(colSums(Y)/sum(Y), 3)
by(Y, seg_vec, function(x) round(colMeans(x), 2))
by(Y, seg_vec, function(x) summary(x))
by(Y, seg_vec, function(x) round(colSums(x)/sum(x), 2))


####EM�A���S���Y���Ő��ݐ��ڍ����������z���f���𐄒�####
##�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z���v�Z���邽�߂̊֐�
LLobz <- function(theta, n, Y, r1, r2, seg){
  LLind <- matrix(0, nrow=nrow(Y), ncol=seg)
  for(i in 1:seg){
    Li <- apply(cbind(n, Y), 1, function(x) dmultinom(x[-1], x[1], theta[i, ]))   #�������z�̖ޓx
    LLind[, i] <- Li
  }
  
  #�ޓx�����������Ă�����A�����Ȑ��𑫂�
  LLind <- ifelse(matrix(apply(LLind, 1, min), nrow=nrow(Y), ncol=seg)==0, LLind+10^-305, LLind)
  
  LLho <- r1 * LLind   #�ϑ��f�[�^�̖ޓx
  z <- LLho / matrix(apply(LLho, 1, sum), nrow=nrow(Y), ncol=seg)   #���ݕϐ�z�̌v�Z
  LLosum <- sum(log(apply(r1 * LLind, 1, sum)))   #�ϑ��f�[�^�̑ΐ��ޓx�̌v�Z
  rval <- list(LLob=LLosum, z=z, LL=LLind)
  return(rval)
}

##�C���f�b�N�X���쐬
index_time <- list()
for(j in 1:seg){
  index_time[[j]] <- subset(1:nrow(ID), ID$time==j)
}


##�����l�̐ݒ�
#�m���̏����l
theta <- matrix(0, nrow=seg, ncol=k)
for(j in 1:seg){
  p0 <- colSums(Y)/sum(Y) + runif(k, 0, 0.75)
  theta[j, ] <- p0 / sum(p0)
}

#�������̏����l
r <- list()
r[[1]] <- c(0.15, 0.2, 0.2, 0.2, 0.25)
R <- r[[1]]

##�ΐ��ޓx�̏�����
y_seg <- Y[index_time[[j]], ]
freq_seg <- freq[index_time[[j]]]
r1_m <- matrix(r[[1]], nrow=n, ncol=seg, byrow=T)

#���ݕϐ�z�Ƒΐ��ޓx�̌v�Z
Z <- matrix(0, nrow=N, ncol=seg)
Z_list <- list()
L <- LLobz(theta, freq_seg, y_seg, r1_m, r1_m, seg)
z <- L$z
Z[index_time[[1]], ] <- z
LL <- L$LLob
LLsum <- LL*seg*2

##���S�f�[�^�̑ΐ��ޓx����p�����[�^���ő剻����
for(j in 1:seg){
  thetaseg <- colSums(matrix(z[, j], nrow=nrow(y_seg), ncol=k) * y_seg) / sum(z[, j] * freq_seg)   #�d�ݕt���������z�̍Ŗސ���
  theta[j, ] <-thetaseg
}

#EM�A���S���Y���̍X�V�X�e�[�^�X
dl <- 100   #EM�X�e�b�v�ł̑ΐ��ޓx�̍��̏����l
tol <- 0.1 
iter <- 0


####EM�A���S���Y���őΐ��ޓx���ő剻####
while(abs(dl) >= tol){   #dl��tol�ȏ�̏ꍇ�͌J��Ԃ�

  #E�X�e�b�v�Ŋϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
  LL <- c()
  for(j in 1:seg){
    y_seg <- Y[index_time[[j]], ]
    freq_seg <- freq[index_time[[j]]]
    
    #���Ԃ�1�Ȃ�1��O�̔����̍����������Ԃ�2�`�Ȃ�1���O�̍����������O���z�Ƃ���
    if(j==1){
      r1 <- r[[j]]
      r1_m <- matrix(r1, nrow=n, ncol=seg, byrow=T)
    } else {
      r1 <- r[[j-1]]
      r1_m <- matrix(r1, nrow=n, ncol=seg, byrow=T)
    }
    
    #�ϑ��f�[�^�̑ΐ��ޓx�Ɛ��ݕϐ�z�̌v�Z
    L <- LLobz(theta, freq_seg, y_seg, r1_m, r1_m, seg)
    Z[index_time[[j]], ] <- L$z
    Z_list[[j]] <- L$z
    LL <- c(LL, L$LLob)
    
    #���������v�Z
    if(j > 1){
      r[[j]] <- colSums(L$z)/n   #�������̌v�Z
    }
  }
  
  #�ΐ��ޓx�̘a�ƑS�̂̍��������v�Z
  LLsum1 <- sum(LL)   
  r[[1]] <- colSums(Z[index_time[[1]], ])/n   #1���ڂ̍�����
  
  ##M�X�e�b�v�Ŋ��S�f�[�^�̑ΐ��ޓx����p�����[�^�𐄒�
  for(j in 1:seg){
    thetaseg <- colSums(matrix(Z[, j], nrow=nrow(Y), ncol=k) * Y) / sum(Z[, j] * freq)   #�d�ݕt���������z�̍Ŗސ���
    theta[j, ] <-thetaseg
  }
  
  ##EM�A���S���Y���̃p�����[�^���X�V
  iter <- iter+1   
  dl <- LLsum-LLsum1
  LLsum <- LLsum1
  print(LLsum1)
}

####���茋�ʂƓK���x�̊m�F####
##�p�����[�^�̌v�Z
round(r_rate <- do.call(rbind, r), 3)   #�������̐���l
z_vec�@<- apply(Z, 1, which.max)   #�Z�O�����g����


##���肳�ꂽ�p�����[�^�Ɛ^�̃p�����[�^�̔�r
round(rbind(theta, P), 3)   #�������z�̃p�����[�^
round(cbind(z_vec, seg_vec))   #�Z�O�����g����
round(cbind(r_rate, r_rate0), 3)   #������


##�K���x�̌v�Z
round(LLsum1, 3)   #�ő剻���ꂽ�ΐ��ޓx
round(AIC <- -2*LLsum1 + 2*(length(theta) + seg), 3)   #AIC
round(BIC <- -2*LLsum1 + log(N)*(length(theta) + seg), 3)   #BIC