######�K�E�X�ߒ���A���f��#####
options(warn=0)
library(MASS)
library(kernlab)
library(GPfit)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)

#set.seed(5698)

####�f�[�^�̔���####
#�f�[�^�̐ݒ�
d <- 1000   #�w�K�f�[�^
k <- 30   #���͕ϐ���
w <- rpois(d, rgamma(d, 5.0, 0.7))

#�e�X�g�f�[�^�̃C���f�b�N�X
index_learn <- 1:d1
index_test <- (d1+1):d

##�f�[�^�̐���
#���n��v�f�̕ϐ��𐶐�
#�T�����𐶐�
week0 <- matrix(diag(7), nrow=d, ncol=7, byrow=T)
week <- week0[, -1]

#�������𐶐�
m_days <- round(runif(ceiling(d/30), 30, 31))
month_vec <- rep(rep(1:12, round(length(m_days)/12))[1:length(m_days)], m_days)
month0 <- matrix(0, nrow=length(month_vec), ncol=12)
for(j in 1:12){
  index <- which(month_vec==j)
  month0[index, j] <- 1
}
month <- month0[1:d, -1]

#�G�ߐ����𐶐�
s_days <- rep(trunc(365/4), ceiling(d/(365/4))) + rbinom(ceiling(d/(365/4)), 1, 1/4)
season_vec <- rep(rep(1:4, length(s_days))[1:length(s_days)], s_days)
season0 <- matrix(0, nrow=length(season_vec), ncol=4)
for(j in 1:4){
  index <- which(season_vec==j)
  season0[index, j] <- 1
}
season <- season0[1:d, -1]


#�������z������͕ϐ��𐶐�
for(rp in 1:1000){
  Data0 <- matrix(0, nrow=d, ncol=k)
  for(j in 1:4){
    index <- which(season0[1:d, j]==1)
    alpha0 <- rep(0.3, k)
    theta <- extraDistr::rdirichlet(1, alpha0)   #�������z�̃p�����[�^
    Data0[index, ] <- rmnom(length(index), w[index], theta)
  }
  if(min(colSums(Data0)) >= 10) break
}

#�f�[�^������
Data <- cbind(week, month, season, Data0)
storage.mode(Data) <- "integer"
sparse_data <- as(Data, "CsparseMatrix")

#�J�[�l���֐��̐���
K <- Data %*% t(Data)

#�K�E�X�ߒ�����Z�O�����g���Ƃɉ����ϐ��𐶐�
sigma <- 2.0
y <- mvrnorm(1, rep(0, d), K + diag(sigma, d))   #�K�E�X�ߒ����牞���ϐ��𐶐�
plot(1:d, y, type="l", xlab="time", ylab="y")


####EM�A���S���Y���ŃK�E�X�ߒ���A���f���𐄒�####
##EM�A���S���Y���̐ݒ�
iter <- 0
LL1 <- -10^100   #�ΐ��ޓx�̏����l
dl <- 100
tol <- 1

#�����l�̐ݒ�
tau1 <- 1
tau2 <- 0.5

#�f�[�^�̐ݒ�
KK <- K %*% K


##EM�A���S���Y���Ńp�����[�^���X�V
while(abs(dl) >= tol){   #dl��tol�ȏ�Ȃ�J��Ԃ�
  #E�X�e�b�v�ŉ�A�W���𐄒�
  beta <- solve(KK + diag(tau1/tau2, d)) %*% t(K) %*% y
  delta <- diag(tau1, d) + tau2*KK
  
  #M�X�e�b�v�Ńn�C�p�[�p�����[�^�𐄒�
  tau1_inv <- (sum(abs(beta^2)) + sum(diag(solve(delta)))) / d
  tau1 <- 1 / tau1_inv
  tau2_inv <- (sum((y - K %*% beta)^2) + sum(diag(KK %*% solve(delta)))) / d
  tau2 <- 1 / tau2_inv
  
  #���Ӗޓx�̍X�V
  diag_tau2 <- diag(tau2_inv, d)
  tau1_KK <- tau1_inv * KK
  Lm <- -1/2*abs(diag_tau2 + tau1_KK) - 1/2*as.numeric((t(y) %*% solve(diag_tau2 + tau1_KK) %*% y))
  LL <- sum(Lm)
  
  ##EM�A���S���Y���̃p�����[�^�̍X�V
  iter <- iter+1
  dl <- LL - LL1
  LL1 <- LL
  print(LL)
}


####���茋�ʂ̊m�F�Ɨv��####
##�\�����ʂƎ����l�̔�r
y_pred <- K %*% beta
plot(1:d, y, type="l", xlab="time", ylab="y", xlim=c(0, d), ylim=c(min(y), max(y)), main="�����l�Ɨ\���l�̎��n��")
par(new=T)
plot(1:d, y_pred, type="l", col=2, xlab="", ylab="", xlim=c(0, d), ylim=c(min(y), max(y)))

##�K���x���m�F
sum((y - y_pred)^2)   #���덷
sum(dnorm(y, y_pred, sd(y - y_pred), log=TRUE))   #�ΐ��ޓx