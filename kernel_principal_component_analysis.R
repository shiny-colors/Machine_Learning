#####�听������#####
library(MASS)
library(plyr)
library(reshape2)
library(kernlab)

####�f�[�^�̔���####
#set.seed(4238)
#���ϗʐ��K���z����̗�������
val <- 10   #�����ϐ��̐�
n <- 1000   #�T���v����

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
corM <- corrM(col=10, lower=0.3, upper=0.9)
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
Sigma <- covmatrix(col=10, corM=corM, lower=20, upper=30)

#�Q���Ƃ̕ϐ��̕��ς��쐬
mu <- c(rnorm(10, 20, 15))

##���ϗʐ��K���z����̗����𔭐�������
val; n
X <- mvrnorm(n=n, mu, Sigma$covariance)

##�f�[�^��v��
round(colMeans(X), 2)   #�ϐ����Ƃ̕���
round(var(X), 2)   #���U�����U�s��
round(cor(X), 2)   #���֍s��
summary(X)   #�v��

##�U�z�}�s��̍쐬
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#�ϐ�1�`5�̎U�z�}�s��
pairs(as.data.frame(X[, 1:5]), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
      upper.panel=panel.cor)
#�ϐ�6�`10�̎U�z�}�s��
pairs(as.data.frame(X[, 6:10]), panel=panel.smooth, bg="lightblue", diag.panel=panel.hist,
      upper.panel=panel.cor)

####�听�����͂����s####
##���f�[�^����听������
Xvar <- var(X)   #���U�����U�s������߂�
PCA <- eigen(Xvar)   #�ŗL�l����(�听������)���s��
PCAval <- PCA$values   #�ŗL�l
PCAvec <- PCA$vectors   #�ŗL�x�N�g��
round(PCAvec %*% diag(PCAval) %*% t(as.matrix(PCAvec)) - Xvar, 3)   #���ْl�������o���Ă��邩����

#PCA�̗v��
round(PCAval/sum(PCAval), 3)   #�������Ƃ̊�^��
round(cumsum(PCAval)/sum(PCAval), 3)   #�ݐϊ�^��

#��3�听���܂Ŏ听�����_���v�Z
PC1 <- X %*% PCAvec[, 1]
PC2 <- X %*% PCAvec[, 2]
PC3 <- X %*% PCAvec[, 3]
PC <- data.frame(PC1, PC2, PC3)
plot(PC, col=4)

####�J�[�l���听������####
##�O�����s��̍쐬
#�������J�[�l��
Xscale <- scale(X)
gram1 <- (1+as.matrix(X) %*% t(as.matrix(X))) + (1+as.matrix(X) %*% t(as.matrix(X)))^2
round(gram1[1:15, 1:15], 3)
round(gram1[985:1000, 985:1000], 3)
round(eigen(gram1)$value, 3)   #������l���ǂ����m�F

#�K�E�X�J�[�l��
sigma <- 1/2
kf_gauss <- function(x1, sigma){
  x1 <- as.matrix(x1)
  g1 <- matrix(t(x1), nrow(x1)^2, ncol(x1), byrow=T)
  g2 <- matrix(rep(x1, c(rep(nrow(x1), nrow(x1)*ncol(x1)))), nrow(x1)^2, ncol(x1))
  g3 <- (g2 - g1)^2
  gram <- exp(-sigma*sqrt(matrix(apply(g3, 1, sum), nrow(x1), nrow(x1), byrow=T)))
  return(gram)
}
gram2 <- kf_gauss(x1=X, sigma=sigma)
round(gram2[1:5, 1:5], 3)
round(gram2[985:1000, 985:1000], 3)
round(eigen(gram2)$value, 3)   #������l���ǂ����m�F

#�֐����g��
L <- kernelMatrix(rbfdot(sigma=1/2), X)   #�K�E�X�J�[�l���ŕϊ�

##�O�����s�񂩂�听�����͂����s
PCA <- eigen(gram1)   #�ŗL�l����(�听������)���s��
PCAval <- PCA$values   #�ŗL�l
PCAvec <- PCA$vectors   #�ŗL�x�N�g��
round(PCAvec %*% diag(PCAval) %*% t(as.matrix(PCAvec)) - gram1, 3)   #���ْl�������o���Ă��邩����

#PCA�̗v��
round(PCAval/sum(PCAval), 3)   #�������Ƃ̊�^��
round(cumsum(PCAval)/sum(PCAval), 3)   #�ݐϊ�^��

#��3�听���܂Ŏ听�����_���v�Z
PC1 <- gram1 %*% PCAvec[, 1]
PC2 <- gram1 %*% PCAvec[, 2]
PC3 <- gram1 %*% PCAvec[, 3]
PC <- data.frame(PC1, PC2, PC3)
plot(PC, col=4)
 
##�֐����g��
kp <- kpca(X, kernel="polydot", kpar=list(degree=2, scale=1))
plot(data.frame(pcv(kp))[, 1:3])