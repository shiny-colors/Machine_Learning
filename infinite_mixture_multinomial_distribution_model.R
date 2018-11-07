#####���������������z���f��#####
library(MASS)
library(Matrix)
library(mclust)
library(reshape2)
library(bayesm)
detach("package:bayesm", unload=TRUE)
library(ExtDist)
library(extraDistr)
library(matrixStats)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(lattice)

#set.seed(1853)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
N <- 10000   #�T���v����
seg <- 10   #�Z�O�����g��
k <- 200   #�ϐ���
w <- rpois(N, rgamma(N, 40, 0.8))   #�T���v��������̕p�x
hist(w, col="grey", breaks=25, main="�A�C�e���w�����̕��z", xlab="�A�C�e���w����")


#�Z�O�����g�̐ݒ�
r <- as.numeric(extraDistr::rdirichlet(1, rep(5.0, seg)))
Z <- rmnom(N, 1, r)
seg_id <- as.numeric(Z %*% 1:seg)


####�����ϐ��̔���####
##�Z�O�����g���ƂɃp�����[�^�̐ݒ�
alpha <- rep(0.2, k)
theta <- thetat <- extraDistr::rdirichlet(seg, alpha)

##�������z���f�[�^�𔭐�
Data <- matrix(0, nrow=N, ncol=k)
for(i in 1:seg){
  index <- which(seg_id==i)
  Data[index, ] <- rmnom(length(index), w[index], theta[i, ])
}
colnames(Data) <- 1:k
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")

####�}���R�t�A�������e�J�����@�Ŗ������������������z���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 10000
burnin <- 1000
keep <- 2
disp <- 20
sbeta <- 1.5
iter <- 0

##���O���z�̐ݒ�
tau <- 1   #�f�B�N�������z�̎��O���z
alpha <- 1   #CRP�̎��O���z
beta <- 1/k

##�����l�̐ݒ�
seg0 <- 2   #�����Z�O�����g��2��
r <- c(0.5, 0.5)   #�������̏����l
par <- rep(0.1, k)   #CRP�p�̃p�����[�^

#�����Z�O�����g��ݒ�
z <- matrix(0, nrow=N, ncol=seg0)
out <- kmeans(Data, seg0)   #kmeans�@

#�Z�O�����g����
z0 <- out$cluster
for(i in 1:seg0){z[z0==i, i] <- 1}   

#�Z�O�����g���Ƃ̃p�����[�^
oldpar0 <- extraDistr::rdirichlet(seg0, par)
oldpar <- (oldpar0+beta) / matrix(rowSums(oldpar0+beta), nrow=seg0, ncol=k, byrow=T)


#�p�����[�^�̊i�[�p�z��
max_seg <- 20
Zi <- matrix(0, nrow=N, ncol=max_seg)
THETA <- array(0, dim=c(max_seg, k, R/keep))
storage.mode(Z) <- "integer"

#�f�[�^�̐ݒ�
const <- lfactorial(rowSums(Data)) - rowSums(lfactorial(Data))


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##�������z�̍����ޓx���v�Z
  #�p�����[�^���Ƃɑΐ��ޓx���v�Z
  LLind0 <- const + Data %*% t(log(oldpar))
  
  #�V�������ݕϐ��̖ޓx�̌v�Z�Ɩޓx�̌���
  par_mean0 <- as.numeric(extraDistr::rdirichlet(1, Data[sample(1:N, 1), ] + beta))
  par_mean <- (par_mean0+beta) / sum(par_mean0+beta)
  LL_new <- const + Data %*% log(par_mean)   #�V�������ݕϐ��̑ΐ��ޓx
  
  LLi0 <- cbind(LLind0, LL_new)
  LLi <- exp(LLi0 - rowMaxs(LLi0))   #�ޓx�ɕϊ�

  ##CRP�̌v�Z
  gamma0 <- cbind(matrix(colSums(z), nrow=N, ncol=ncol(z), byrow=T) - z, alpha)
  gamma1 <- LLi * gamma0/(N-1-alpha)
  
  ##�������z�����ݕϐ����T���v�����O
  z_rate <- gamma1 / rowSums(gamma1)   #���ݕϐ�z�̊����m��
  z <- rmnom(N, 1, z_rate)
  z <- z[, colSums(z) > 0]
  
  ##�������z�̃p�����[�^���X�V
  dir_par <- t(t(Data) %*% z) + tau   #�f�B�N�������z�̃p�����[�^
  oldpar <- extraDistr::rdirichlet(ncol(z), dir_par)   #�������z����p�����[�^�𐶐�

  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  if(rp%%keep==0){
    mkeep <- rp/keep
    THETA[1:nrow(oldpar), , mkeep] <- oldpar
    
    #�J��Ԃ������o�[���C���𒴂�����p�����[�^���i�[
    if(rp >= burnin){
      Zi[, 1:ncol(z)] <- Zi[, 1:ncol(z)] + z
    }   
    
    if(rp%%disp==0){
      print(rp)
      print(colSums(z))
    }
  }
}

####�T���v�����O���ʂ̉����Ɨv��####
#�o�[���C������
burnin1 <- R/(keep+2)   
burnin2 <- 1000

##�T���v�����O���ʂ��v���b�g
matplot(t(THETA[1, , 1:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(THETA[2, , 1:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(THETA[3, , 1:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(THETA[4, , 1:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(THETA[5, , 1:(R/keep)]), type="l", ylab="�p�����[�^")
matplot(t(THETA[6, , 1:(R/keep)]), type="l", ylab="�p�����[�^")

##�T���v�����O���ʂ̎��㕽��
mcmc_seg <- sum(colSums(Zi) >= N)   #���肳�ꂽ�Z�O�����g��

#���ݕϐ�z�̐����
round(Z_mu <- (Zi/rowSums(Zi))[, colSums(Zi) > 0], 3)   #���ݕϐ��̊����m��
colnames(Z_mu) <- 1:ncol(Z_mu)
round(colMeans(Z_mu), 3)   #������

#�������z�̃p�����[�^�̐����
theta_mu <- matrix(0, nrow=mcmc_seg, ncol=k)
for(i in 1:mcmc_seg){
  theta_mu[i, ] <- colMeans(t(THETA[i, , burnin1:(R/keep)]))
}
round(cbind(t(theta_mu), t(thetat)), 3) #�^�̃p�����[�^�Ɣ�r
