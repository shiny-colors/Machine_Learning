#####���������K�E�X���z���f��#####
library(MASS)
library(mclust)
library(reshape2)
library(gtools)
library(bayesm)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

####���ϗʐ��K���z�̗����𔭐�������֐����`####
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper, eigen_lower, eigen_upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, runif(1, eigen_lower, eigen_upper), Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}


##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
  }
  diag(c) <- m
  cc <- c * corM
  #�ŗL�l�����ŋ����I�ɐ���l�s��ɏC������
  UDU <- eigen(cc)
  val <- UDU$values
  vec <- UDU$vectors
  D <- ifelse(val < 0, val + abs(val) + 0.00001, val)
  covM <- vec %*% diag(D) %*% t(vec)
  data <- list(covM, cc,  m)
  names(data) <- c("covariance", "cc", "mu")
  return(data)
}

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
hh <- 5000
seg <- 5
k <- 4

#�Z�O�����g�����̐ݒ�
seg_id <- rep(1:seg, rep(hh/seg, seg))

##�p�����[�^�̐ݒ�
#���ύ\���̔���
mu0 <- matrix(runif(k*seg, -5, 5), nrow=seg, ncol=k)

#���U�����U�s��̔���
Cov0 <- list()
for(j in 1:seg){
  Cor <- corrM(k, -0.55, 0.7, 0.1, 0.3)
  Cov0[[j]] <- covmatrix(k, Cor, 0.5, 2)$covariance
}

##���ϗʐ��K���z����f�[�^�𔭐�
Data <- matrix(0, nrow=hh, ncol=k)
for(j in 1:seg){
  Data[seg_id==j, ] <- mvrnorm(length(seg_id[seg_id==j]), mu0[j, ], Cov0[[j]])
}

#��ϗʂ��Ƃ̌��ʂ�����
plot(Data[, 1:2], col=seg_id[seg_id %in% 1:5], pch=20, xlab="�f�[�^1�̒l", ylab="�f�[�^2�̒l", main="������ϗʐ��K���z�̃v���b�g")
plot(Data[, c(1, 3)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="�f�[�^1�̒l", ylab="�f�[�^3�̒l", main="������ϗʐ��K���z�̃v���b�g")
plot(Data[, c(1, 4)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="�f�[�^1�̒l", ylab="�f�[�^4�̒l", main="������ϗʐ��K���z�̃v���b�g")
plot(Data[, 2:3], col=seg_id[seg_id %in% 1:5], pch=20, xlab="�f�[�^2�̒l", ylab="�f�[�^3�̒l", main="������ϗʐ��K���z�̃v���b�g")
plot(Data[, c(2, 4)], col=seg_id[seg_id %in% 1:5], pch=20, xlab="�f�[�^2�̒l", ylab="�f�[�^4�̒l", main="������ϗʐ��K���z�̃v���b�g")
plot(Data[, 3:4], col=seg_id[seg_id %in% 1:5], pch=20, xlab="�f�[�^3�̒l", ylab="�f�[�^4�̒l", main="������ϗʐ��K���z�̃v���b�g")


####�}���R�t�A�������e�J�����@�Ŗ��������K�E�X���z���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
keep <- 2
sbeta <- 1.5
iter <- 0

##���ϗʐ��K���z�̖ޓx�֐�
dmv <- function(x, mean.vec, S, S_det, S_inv){
  LLo <- (2*pi)^(-nrow(S)/2) * S_det^(-1/2) *
    exp(-1/2 * (x - mean.vec) %*% S_inv %*% (x - mean.vec))
  return(LLo)
}

##���O���z�̐ݒ�
#���ϗʐ��K���z�̎��O���z
mean0 <- rep(0, k)
sigma0 <- diag(100, k)
sigma0_inv <- solve(sigma0)
mean_z <- colMeans(Data)   #��ĕ��z�̕���
sigma_z <- diag(diag(var(Data))) * 7.5   #��ĕ��z�̕��U

#�t�E�B�V���[�g���z�̎��O���z
nu <- k+1
V <- nu * diag(k)

#�f�B�N�������z�̎��O���z
alpha <- 1

##�����l�̐ݒ�
seg0 <- 2   #�����Z�O�����g��2��
res <- Mclust(Data, 2)   #�������K���z�𐄒�

#���ݕϐ�z�̏����l
z_vec <- apply(res$z, 1, which.max)
z <- matrix(0, nrow=hh, ncol=seg0)
for(i in 1:hh) {z[i, z_vec[i]] <- 1}

#��A�p�����[�^�̏����l
oldmean <- t(res$parameters$mean)
oldcov <- res$parameters$variance$sigma

##�p�����[�^�̊i�[�p�z��
max_seg <- 15
Z <- matrix(0, nrow=R/keep, ncol=max_seg)
Mu <- array(0, dim=c(max_seg, k, R/keep))
Cov <- list()
storage.mode(Z) <- "integer"


####MCMC�Ńp�����[�^���T���v�����O
for(rp in 1:R){
  
  ##���ݕϐ�z���T���v�����O
  z_len <- length(unique(z_vec))
  LLi <- matrix(0, nrow=hh, ncol=z_len+1)
  
  #�T���v�����O�ς݂̐��ݕϐ��̖ޓx
  for(j in 1:z_len){
    LLi[, j] <- dmvnorm(Data, oldmean[j, ], oldcov[, , j])
  }
  #�V�������ݕϐ��̖ޓx
  LLi[, z_len+1] <- dmvnorm(Data, mean_z, sigma_z)
  
  #CRP���v�Z
  gamma0 <- cbind(matrix(colSums(z), nrow=hh, ncol=z_len, byrow=T) - z, alpha)
  gamma1 <- LLi * gamma0/(hh-1-alpha)
  
  #�������z�����ݕϐ�z���T���v�����O
  z <- t(apply(gamma1/rowSums(gamma1), 1, function(x) rmultinom(1, 1, x)))
  z <- z[, colSums(z) > 0]
  z_vec <- z %*% 1:ncol(z)
  
  ##���ϗʐ��K���z�̃p�����[�^���T���v�����O
  #��������Ȃ�����z�͏������Ă���
  z_cnt <- length(colSums(z) > 0)
  if(z_cnt > nrow(oldmean)){
    oldmean <- rbind(oldmean, 0)
  }
  oldmean <- oldmean[1:z_cnt, ]
  oldcov <- array(0, dim=c(k, k, z_cnt))
  
  #���ϗʐ��K���z�̕��ς���ѕ��U�����U�s����M�u�X�T���v�����O
  for(j in 1:z_cnt){
    
    #�C���f�b�N�X���쐬
    index <- subset(1:nrow(z), z[, j]==1)
    
    #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
    Vn <- nu + length(index)
    er <- Data[index, ] - matrix(oldmean[j, ], nrow=length(index), ncol=ncol(Data), byrow=T)
    R_par <- solve(V) + t(er) %*% er
    
    oldcov[, , j] <- rwishart(Vn, solve(R_par))$IW   #�t�E�B�V���[�g���z���番�U�����U�s����T���v�����O
    
    #���ς��T���v�����O
    if(length(index) > 1){
      mean_mu <- length(index)/(1+length(index))*colMeans(Data[index, ])   #���σp�����[�^�̕���
    } else {
      mean_mu <- length(index)/(1+length(index))*Data[index, ] 
    }
    mean_cov <- oldcov[, , j] / (1+length(index))   #���σp�����[�^�̕��U�����U�s��
    oldmean[j, ] <- mvrnorm(1, mean_mu, mean_cov)
  }
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    Mu[1:nrow(oldmean), , mkeep] <- oldmean
    Cov[[mkeep]] <- oldcov
    if(rp >= R/2){Z[, 1:ncol(z)] <- Z[, 1:ncol(z)] + z}   #�J��Ԃ������ő唽�����̔����𒴂�����p�����[�^���i�[
    
    print(rp)
    print(colSums(z))
  }
}

####�T���v�����O���ʂ̗v��Ɖ���####
#�o�[���C������
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##�T���v�����O���ʂ��v���b�g
matplot(t(Mu[1, , burnin2:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(Mu[2, , burnin2:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(Mu[3, , burnin2:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(Mu[4, , burnin2:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(Mu[5, , burnin2:(R/keep)]), type="l", ylab="�p�����[�^")

##�T���v�����O���ʂ̎��㕽��
mcmc_seg <- sum(colSums(Z) > 10000)   #���肳�ꂽ�Z�O�����g��

#���ݕϐ�z�̐����
round(Z_mu <- (Z/rowSums(Z))[, colSums(Z) > 0], 3)   #���ݕϐ��̊����m��
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #������

#���ς̐����
mean_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  mean_mu[i, ] <- colMeans(t(Mu[i, , burnin1:(R/keep)]))
}
round(rbind(mean_mu, mu0), 3)   #�^�̃p�����[�^�Ɣ�r

#���U�����U�s��̐����
cov_mu0 <- Cov[[burnin1]]
for(i in (burnin1+1):(R/keep)){
  cov_mu0 <- cov_mu0 + Cov[[i]][, , 1:mcmc_seg]
}
round(Cov_mu <- cov_mu0/length(burnin1:(R/keep)), 3)   #���U�����U�s��̎��㕽��
lapply(Cov0, function(x) round(x, 3))   #�^�̕��U�����U�s��
