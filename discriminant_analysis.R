#####���ʃ��f��#####
library(MASS)
library(plyr)
library(reshape2)
####���Q�̐������ʃ��f��####
####�f�[�^�̔���####
#set.seed(4235)
#���ϗʐ��K���z����̗�������
k <- 4   #�Q�̐�
val <- 6   #�����ϐ��̐�
n <- 400   #�Q���Ƃ̊w�K�Ɏg�����߂̃f�[�^��
nt <- 300   #�Q���Ƃ̃e�X�g�Ɏg�����߂̃f�[�^��

##���ϗʐ��K���z����̗����𔭐�������
#�C�ӂ̑��֍s������֐����`
corrM <- function(col, lower, upper){
  diag(1, col, col)
  
  rho <- matrix(runif(col^2, lower, upper), col, col)
  rho[upper.tri(rho)] <- 0
  Sigma <- rho + t(rho)
  diag(Sigma) <- 1
  Sigma
  (X.Sigma <- eigen(Sigma))
  (Lambda <- diag(X.Sigma$values))
  P <- X.Sigma$vector
  P %*% Lambda %*% t(P)
  
  #�V�������֍s��̒�`�ƑΊp������1�ɂ���
  (Lambda.modified <- ifelse(Lambda < 0, 10e-6, Lambda))
  x.modified <- P %*% Lambda.modified %*% t(P)
  normalization.factor <- matrix(diag(x.modified),nrow = nrow(x.modified),ncol=1)^0.5
  Sigma <- x.modified <- x.modified / (normalization.factor %*% t(normalization.factor))
  eigen(x.modified)
  diag(Sigma) <- 1
  round(Sigma, digits=3)
  return(Sigma)
}

#�Q���Ƃ̑��֍s����쐬(�Q�ł��ׂē���)
corM <- corrM(col=6, lower=-0.2, upper=0.2)
eigen(corM)

##���֍s�񂩂番�U�����U�s����쐬����֐����`
covmatrix <- function(col, corM, lower, upper){
  m <- abs(runif(col, lower, upper))
  c <- matrix(0, col, col)
  for(i in 1:col){
    for(j in 1:col){
      c[i, j] <- sqrt(m[i]) * sqrt(m[j])
    }
    diag(c) <- m   #�Ίp�s������̕��U�ɖ߂�
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

#���U�����U�s����쐬(�Q�ł��ׂē���)
Sigma1 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma2 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma3 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
Sigma4 <- covmatrix(col=6, corM=corM, lower=10, upper=20)
S <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance, Sigma4$covariance)

#�Q���Ƃ̕ϐ��̕��ς��쐬
mu1 <- c(rnorm(6, 12, 10))
mu2 <- c(rnorm(6, 18, 10))
mu3 <- c(rnorm(6, 6, 10))
mu4 <- c(rnorm(6, 24, 10))
mu <- list(mu1, mu2, mu3, mu4)

##���ϗʐ��K���z����̗����𔭐�������
k; n; nt
X <- matrix(0, 0, 6)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], S[[kk]])
  X <- rbind(X, xx)
}

#���t�f�[�^�̃x�N�g�����쐬
y <- rep(1:4, rep(700, 4))

#�f�[�^������
YX <- as.data.frame(cbind(y, X))
by(YX, YX[, 1] , function(x) summary(x))   #�Q���Ƃ̗v��֐�
by(YX[, 2:7], YX[, 1] , function(x) cor(x))   #�Q���Ƃ̑���
plot(YX[, 2:7], col=YX[, 1])   #�U�z�}

#�e�X�g�f�[�^�Ɗw�K�f�[�^�𕪗�
#�w�K�f�[�^
YXl <- rbind(YX[1:400, ], YX[701:1100, ], YX[1401:1800, ], YX[2101:2500, ])
table(YXl[, 1])

#�e�X�g�f�[�^
YXt <- rbind(YX[401:700, ], YX[1101:1400, ], YX[1801:2100, ], YX[2501:2800, ])
table(YXt[, 1])

####�������ʃ��f���Ŋw�K�f�[�^���番�ފ���w�K����####
#�Q������
mucat <- matrix(0, k, val)
for(i in 1:k){
  muc <- apply(YXl[YXl[, 1]==i, 2:7], 2, mean)
  mucat[i, ] <- (muc)
}

#�Q�S�̂̕���
muall <- colMeans(YXl[, 2:7])
names(muall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")

#�Q�S�̂̕��U�����U�s��
covall <- var(YXl[, 2:7])
rownames(covall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")
colnames(covall)[1:6] <- c("v1", "v2", "v3", "v4", "v5", "v6")

#�Q�Ԃ̕��U�����U�s��
covB <- 1/k * t(mucat - muall) %*% (mucat - muall)

##�ŗL�l���������ČQ�Ԃ̕����x���ő剻������𓾂�
covA <- solve(covall)
M <- eigen(covA %*% as.matrix(covB))   #�ŗL�l��������
R <- M$values   #�s��̃����N�͌Q��-1
a <- M$vectors   #�ŗL�x�N�g�������ʊ֐��̌W��
a
(cont <- R/sum(R))   #��^��
(cumcont <-cumsum(cont))   #�ݐϊ�^��

#��2�ŗL�l�܂ł�p���āA2������ԏ�Ɏˉe�����s����v���b�g
#�w�K�f�[�^�ɑ΂��鍇���ϗʂƃv���b�g
Zl <- as.matrix(YXl[, 2:7]) %*% t(a[1:3, ])
plot(as.data.frame(Zl), col=YXl[, 1])

#�e�X�g�f�[�^�ɑ΂��鍇���ϗʂƃv���b�g
Zt <- as.matrix(YXt[, 2:7]) %*% t(a[1:3, ])
plot(as.data.frame(Zt))

####���ʐ��\���m�F####
#�w�K�f�[�^�̐��������v�Z
zm <- matrix(0, k, 3)
for(i in 1:k){
  zm[i, ] <- apply(Zl[YXl[, 1]==i, ], 2, mean)
}

#���ʌ��ʂ����߂�
D <- matrix(0, nrow=nrow(ZZl), ncol=k)
for(i in 1:k){
  ZZl <- Zl - matrix(zm[i, ], nrow=nrow(Zl), ncol=ncol(Zl), byrow=T)
  d <- matrix(apply(ZZl, 1, function(x) t(x) %*% x), nrow=nrow(ZZl), ncol=1)
  D[, i] <- d
}
res <- apply(D, 1, which.min)   #���ʌ���
pre <- data.frame(correct=YXl[, 1], predict=res)   #�����f�[�^�ƌ���
table(pre)   #�딻�ʕ\
round(apply(table(pre), 1, function(x) x/sum(x)), 3)   #�Q���Ƃ̐�����
sum(diag(table(pre)))/sum(table(pre))   #������


####2������####
##�قȂ鋤���U�s��̑��ϗʐ��K���z���痐������
#�Q���Ƃ̑��֍s����쐬(�Q�ňقȂ�)
k <- 3   #�Q��
val <- 6   #�����ϐ���
n <- 500   #�Q���Ƃ̊w�K�f�[�^
nt <- 500   #�Q���Ƃ̃e�X�g�f�[�^
N <- 3000   #�S�T���v����

corM1 <- corrM(col=6, lower=-0.25, upper=0.3)
corM2 <- corrM(col=6, lower=-0.3, upper=0.4)
corM3 <- corrM(col=6, lower=-0.4, upper=0.55)

#���U�����U�s����쐬(�Q�ł��ׂē���)
Sigma1 <- covmatrix(col=6, corM=corM1, lower=15, upper=26)
Sigma2 <- covmatrix(col=6, corM=corM2, lower=18, upper=30)
Sigma3 <- covmatrix(col=6, corM=corM3, lower=12, upper=18)
S <- list(Sigma1$covariance, Sigma2$covariance, Sigma3$covariance)

#�Q���Ƃ̕ϐ��̕��ς��쐬
mu1 <- c(rnorm(6, 12, 10))
mu2 <- c(rnorm(6, 18, 10))
mu3 <- c(rnorm(6, 8, 10))
mu <- list(mu1, mu2, mu3)

##���ϗʐ��K���z����̗����𔭐�������
X <- matrix(0, 0, 6)
for(kk in 1:k){
  xx <- mvrnorm(n=n+nt, mu[[kk]], S[[kk]])
  X <- rbind(X, xx)
}

#���t�f�[�^�̃x�N�g�����쐬
y <- rep(1:3, rep(1000, 3))

#�f�[�^������
YX <- data.frame(y, X)
by(YX, YX[, 1] , function(x) summary(x))   #�Q���Ƃ̗v��֐�
by(YX[, 2:7], YX[, 1] , function(x) cor(x))   #�Q���Ƃ̑���
plot(YX[, 2:7], col=YX[, 1])   #�U�z�}

#�e�X�g�f�[�^�Ɗw�K�f�[�^�𕪗�
#�w�K�f�[�^
YXl <- rbind(YX[1:500, ], YX[1001:1500, ], YX[2001:2500, ])
table(YXl[, 1])

#�e�X�g�f�[�^
YXt <- rbind(YX[501:1000, ], YX[1501:2000, ], YX[2501:3000, ])
table(YXt[, 1])

####2�����ʂ̐���####
##�Q���Ƃ̕��σx�N�g���ƕ��U�����U�s����v�Z
#�Q���Ƃ̕��σx�N�g��
mv <- matrix(0, k, val)
for(mm in 1:k){
  m <- apply(YXl[YXl[, 1]==mm, 2:7], 2, mean)
  mv[mm, ] <- m
}

#�Q���Ƃ̕��U�����U�s��
Sk <- list()
for(ss in 1:k){
  m <- var(YXl[YXl[, 1]==ss, 2:7])
  Sk[[ss]] <- m
}

##�T���v�����Ƃɂ��ꂼ��̌Q�̃}�n���m�r�X�̔ċ��������߂Ă����Ƃ��������Q�ɏ���������
##�w�K�f�[�^�ɂ��ċ��߂�
HL <- matrix(0, nrow=nrow(YXl), ncol=k)
S <- list(solve(Sk[[1]]), solve(Sk[[2]]), solve(Sk[[3]]))
for(kk in 1:k){
  h <- apply(YXl[, 2:7], 1, function(x) t(x-mv[kk, ]) %*% S[[kk]] %*% (x-mv[kk, ]))
  HL[, kk] <- h
}

#�w�K�f�[�^�̔��ʌ���
pred <- apply(HL, 1, which.min)   #���ʌ���
res <- data.frame(true=YXl[, 1], pred)   #�^�̌Q�ƌ���
table(res)   #�딻�ʕ\
sum(diag(table(res)))/sum(table(res))   #�딻�ʗ�


##�e�X�g�f�[�^�ɂ��ċ��߂�
HT <- matrix(0, nrow=nrow(YXt), ncol=k)
S <- list(solve(Sk[[1]]), solve(Sk[[2]]), solve(Sk[[3]]))
for(kk in 1:k){
  h <- apply(YXt[, 2:7], 1, function(x) t(x-mv[kk, ]) %*% S[[kk]] %*% (x-mv[kk, ]))
  HT[, kk] <- h
}

#�e�X�g�f�[�^�̔��ʌ���
pred <- apply(HT, 1, which.min)   #���ʌ���
res <- data.frame(true=YXt[, 1], pred)   #�^�̌Q�ƌ���
table(res)   #�딻�ʕ\
sum(diag(table(res)))/sum(table(res))   #�딻�ʗ�