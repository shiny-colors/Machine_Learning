#####Hierarchical Hidden Marcov Mixture Model#####
options(warn=2)
library(MASS)
library(lda)
library(RMeCab)
library(matrixStats)
library(Matrix)
library(bayesm)
library(HMM)
library(extraDistr)
library(reshape2)
library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)

'%!in%' <- function(a,b) ! a %in% b

#set.seed(5723)

####�f�[�^�̔���####
##�f�[�^�̐ݒ�
k1 <- 4   #������
k2 <- 7   #��Ԑ�
k3 <- 10   #�g�s�b�N��
d <- 5000   #������
s <- rep(300, k1)   #�������Ƃ̌�b��
v <- sum(s)   #����b��
v1 <- c(1, (cumsum(s)+1)[-k1])
v2 <- cumsum(s)
w <- rpois(d, rgamma(d, 70, 0.5))   #�������Ƃ̒P�ꐔ
f <- sum(w)   #����b��

##ID�̐ݒ�
d_id <- rep(1:d, w)
t_id <- as.numeric(unlist(tapply(1:f, d_id, rank)))
index_end <- as.numeric(tapply(1:f, d_id, max))

##�p�����[�^�̎��O���z�̐ݒ�
#�}���R�t���ڊm���̃f�B���N�����z�̎��O���z��ݒ�
alpha11 <- rep(5.0, k1)
alpha12 <- matrix(1.5, nrow=k1, ncol=k1)
diag(alpha12) <- 10^-300
alpha21 <- rep(1.5, k2)
alpha22 <- matrix(0.2, nrow=k2, ncol=k2)
diag(alpha22) <- 1.0

#�؂�ւ��m���̃x�[�^���z�̎��O���z��ݒ�
s0 <- 7.5
v0 <- 50.0

#�g�s�b�N���z�̃f�B���N�����z�̎��O���z��ݒ�
alpha3 <- array(0.15, dim=c(k2, k3, k1))

#�g�s�b�N���Ƃ̒P�ꕪ�z�̃f�B���N�����z�̎��O���z��ݒ�
gamma <- matrix(0.0001, nrow=k1, ncol=v)
for(j in 1:k1){
  gamma[j, v1[j]:v2[j]] <- 0.1
}

##���ׂĂ̒P�ꂪ���������܂Ń��[�v
rp <- 0
repeat {
  rp <- rp + 1
  print(rp)
  
  ##���O���z����p�����[�^�𐶐�
  #�f�B���N�����z����}���R�t���ڊm���𐶐�
  theta11 <- thetat11 <- as.numeric(extraDistr::rdirichlet(1, alpha11))
  theta12 <- thetat12 <- extraDistr::rdirichlet(k1, alpha12)
  theta21 <- thetat21 <- extraDistr::rdirichlet(k1, alpha21)
  theta22 <- array(0, dim=c(k2, k2, k1))
  for(j in 1:k1){
    theta22[, , j] <- extraDistr::rdirichlet(k2, alpha22)
  }
  thetat22 <- theta22
  
  #�x�[�^���z����؂�ւ��ϐ��𐶐�
  beta <- betat <- rbeta(d, s0, v0)
  
  #�f�B���N�����z����g�s�b�N���z�𐶐�
  theta3 <- array(0, dim=c(k2, k3, k1))
  for(j in 1:k1){
    theta3[, , j] <- extraDistr::rdirichlet(k2, alpha3[, , j])
    
  }
  thetat3 <- theta3
  
  #�f�B���N�����z����P�ꕪ�z�𐶐�
  phi <- array(0, dim=c(k3, v, k1))
  for(j in 1:k1){
    repeat {
      phi[, , j] <- extraDistr::rdirichlet(k3, gamma[j, ])
      if(min(colMaxs(phi[, v1[j]:v2[j], j])) > (k3^2*2)/f) break
    }
  }
  phit <- phi
  
  ##HHMM���f���̉��肩��f�[�^�𐶐�
  #�f�[�^�̊i�[�p
  z1_list <- list()
  z21_list <- list()
  z22_list <- list()
  z3_list <- list()
  wd_list <- list()
  WX <- matrix(0, nrow=d, ncol=v)
  
  for(i in 1:d){
    
    #�؊����ϐ��𐶐�
    z1_vec <- rbinom(w[i], 1, beta[i])
    z1_vec[1] <- 1
    
    ##�������z���}���R�t��Ԑ��ڂ𐶐�
    z21_vec <- rep(0, w[i])
    z22_vec <- rep(0, w[i])
    z3_vec <- rep(0, w[i])
    word_vec <- rep(0, w[i])
    words <- matrix(0, nrow=w[i], ncol=v)
    
    for(j in 1:w[i]){
      if(j==1){
        #��ʊK�w�̏�Ԑ��ڂ𐶐�
        z21 <- rmnom(1, 1, theta11)
        z21_vec[j] <- as.numeric(z21 %*% 1:k1)
        
        #���ʊK�w�̏�Ԑ��ڂ𐶐�
        z22 <- rmnom(1, 1, theta21[z21_vec[j], ])
        z22_vec[j] <- as.numeric(z22 %*% 1:k2)
        
      } else {
        
        if(z1_vec[j]==1){
          #��ʊK�w�̏�Ԑ��ڂ𐶐�
          z21 <- rmnom(1, 1, theta12[z21_vec[j-1], ])
          z21_vec[j] <- as.numeric(z21 %*% 1:k1)
          
          #���ʊK�w�̏�Ԑ��ڂ𐶐�
          z22 <- rmnom(1, 1, theta21[z21_vec[j], ])
          z22_vec[j] <- as.numeric(z22 %*% 1:k2)
          
        } else {
          
          #��ʊK�w�̏�Ԑ��ڂ𐶐�
          z21_vec[j] <- z21_vec[j-1]
          
          #���ʊK�w�̏�Ԑ��ڂ𐶐�
          z22 <- rmnom(1, 1, theta22[z22_vec[j-1], , z21_vec[j]])
          z22_vec[j] <- as.numeric(z22 %*% 1:k2)
        }
      }
      ##�g�s�b�N�ƒP��𐶐�
      #�������z����g�s�b�N�𐶐�
      z3 <- rmnom(1, 1, theta3[z22_vec[j], , z21_vec[j]])
      z3_vec[j] <- as.numeric(z3 %*% 1:k3)
      
      #�g�s�b�N����P��𐶐�
      word <- rmnom(1, 1, phi[z3_vec[j], , z21_vec[j]])
      words[j, ] <- word
      word_vec[j] <- as.numeric(word %*% 1:v)
    }
    
    #���������f�[�^���i�[
    z1_list[[i]] <- z1_vec
    z21_list[[i]] <- z21_vec
    z22_list[[i]] <- z22_vec
    z3_list[[i]] <- z3_vec
    wd_list[[i]] <- word_vec
    WX[i, ] <- colSums(words)
  }
  if(min(colSums(WX)) > 0){
    break 
  }
}

##�f�[�^��ϊ�
Z1 <- unlist(z1_list)
Z21 <- unlist(z21_list)
Z22 <- unlist(z22_list)
Z3 <- unlist(z3_list)
wd <- unlist(wd_list)
sparse_data <- sparseMatrix(1:f, wd, dims=c(f, v))
sparse_data_T <- t(sparse_data)


####�}���R�t�A�������e�J�����@��HHMM���f���𐄒�####
##�A���S���Y���̐ݒ�
R <- 3000
keep <- 2  
iter <- 0
burnin <- 1000
disp <- 10

##���O���z�̐ݒ�
#�n�C�p�[�p�����[�^�̎��O���z
s0 <- 0.5
v0 <- 0.5
alpha01 <- 1.0
alpha02 <- 1.0
beta01 <- 0.1
beta02 <- 10^-10

##MCMC�p�C���f�b�N�X���쐬
#�P�ꏇ���̃C���f�b�N�X
max_word <- max(t_id)
index_t11 <- which(t_id==1)
index_t21 <- index_t22 <- list()
n_vec <- c(length(index_t11), rep(0, max_word-1))
for(j in 1:max_word){
  index_t21[[j]] <- which(t_id==j)-1
  index_t22[[j]] <- which(t_id==j)
  n_vec[j] <- length(index_t22[[j]])
}
vec_k1 <- rep(1, k1); vec_k2 <- rep(1, k2); vec_k3 <- rep(1, k3)

#��������ђP��̃C���f�b�N�X
doc_list <- doc_vec <- list()
wd_list <- wd_vec <- list()
for(i in 1:d){
  doc_list[[i]] <- which(d_id==i)
  doc_vec[[i]] <- rep(1, length(doc_list[[i]]))
}


##�^�l��ݒ�
#�p�����[�^�̐^�l
const <- lfactorial(w) - rowSums(lfactorial(WX))   #�������z�̖��x�֐��̑ΐ��ޓx�̒萔
beta <- betat
theta11 <- thetat11
theta12 <- thetat12
theta21 <- thetat21
theta22 <- thetat22
theta3 <- thetat3
phi <- phit


##�����l��ݒ�
#�p�����[�^�̏����l�𐶐�
const <- lfactorial(w) - rowSums(lfactorial(WX))   #�������z�̖��x�֐��̑ΐ��ޓx�̒萔
beta <- rep(0.5, d)
theta11 <- rep(1/k1, k1)
theta12 <- matrix(1/(k1-1), nrow=k1, ncol=k1); diag(theta12) <- 0
theta21 <- extraDistr::rdirichlet(k1, rep(10.0, k2))
theta22 <- array(0, dim=c(k2, k2, k1))
theta3 <- array(0, dim=c(k2, k3, k1))
phi <- array(0, dim=c(k3, v, k1))
for(j in 1:k1){
  theta22[, , j] <- extraDistr::rdirichlet(k2, rep(10.0, k2))
  theta3[, , j] <- extraDistr::rdirichlet(k2, rep(10.0, k3))
  phi[, , j] <- extraDistr::rdirichlet(k3, rep(10.0, v))
}

##�p�����[�^�̕ۑ��p�z��
BETA <- matrix(0, nrow=R/keep, ncol=d)
THETA11 <- matrix(0, nrow=R/keep, k1)
THETA12 <- array(0, dim=c(k1, k1, R/keep))
THETA21 <- array(0, dim=c(k1, k2, R/keep))
THETA22 <- array(0, dim=c(k2, k2, k1, R/keep))
PHI <- array(0, dim=c(k3, v, k1, R/keep))
SEG1 <- rep(0, f)
SEG21 <- matrix(0, nrow=f, ncol=k1)
SEG22 <- matrix(0, nrow=f, ncol=k2)
SEG3 <- matrix(0, nrow=f, ncol=k3)

#�ΐ��ޓx�̊�l
LLst <- sum(sparse_data %*% log(colSums(sparse_data) / sum(sparse_data)))
Z21_freq <- as.numeric((plyr::count(Z21))$freq)


####MCMC�Ńp�����[�^���T���v�����O####
for(rp in 1:R){
  
  ##���ݕϐ��̔z��̐ݒ�
  Zi1 <- z21_vec <- z22_vec <- rep(0, f)
  Zi21 <- matrix(0, nrow=f, ncol=k1)
  Zi22 <- matrix(0, nrow=f, ncol=k2)
  
  ##�g�s�b�N���Ƃ̊��Җޓx�𐄒�
  LLt <- array(0, dim=c(f, k2, k1))
  for(j in 1:k1){
    LLt[, , j] <- t(phi[, , j])[wd, ] %*% t(theta3[, , j]) 
  }
  
  ##1�P��ڂ̏�ʊK�w�̐��ݏ�Ԃ𐶐�
  #��ʊK�w�̊��Җޓx
  LLt21 <- LLt[index_t11, , ]
  par <- matrix(0, nrow=d, ncol=k1)
  for(j in 1:k1){
    par[, j] <- LLt21[, , j] %*% theta21[j, ]
  }
  
  #�������z�����ʊK�w�̏�Ԃ𐶐�
  r <- matrix(theta11, nrow=d, ncol=k1, byrow=T)
  par_rate <- r*par / as.numeric(par %*% theta11)   #���ݕϐ��̊����m��
  Zi21[index_t11, ] <- rmnom(d, 1, par_rate)   #�������z�����Ԃ𐶐�
  z21_vec[index_t11] <- as.numeric(Zi21[index_t11, ] %*% 1:k1)
  
  ##1�P��ڂ̉��ʊK�w�̐��ݏ�Ԃ𐶐�
  #��ʊK�w�̏�Ԃɉ����Ėޓx�v�Z
  zi21 <- Zi21[index_t11, ]
  par <- matrix(0, nrow=d, ncol=k2)
  for(j in 1:k1){
    par <- par + LLt21[, , j] * zi21[, j]
  }
  
  #�������z���牺�ʊK�w�̏�Ԃ𐶐�
  r <- theta21[z21_vec[index_t11], ]
  par_rate <- r*par / (rowSums2(r*par))   #���ݕϐ��̊����m��
  Zi22[index_t11, ] <- rmnom(d, 1, par_rate)   #�������z�����Ԃ𐶐�
  z22_vec[index_t11] <- as.numeric(Zi22[index_t11, ] %*% 1:k2)
    
  
  ##2�P��ڈȍ~�̐��ݏ�Ԃ𐶐�
  #���ڍs��̔z���ݒ�
  rf21 <- matrix(0, nrow=k1, ncol=k1); rf22 <- array(0, dim=c(k2, k2, k1))
  for(pd in 2:max_word){
    
    ##���ݏ�Ԃ��؂�ւ�邩�ǂ����𐶐�
    #�f�[�^�̐ݒ�
    index1 <- index_t21[[pd]]; index2 <- index_t22[[pd]]
    z21_vec_j <- z21_vec[index1]; zi21_j <- Zi21[index1, , drop=FALSE]
    z22_vec_j <- z22_vec[index1]; zi22_j <- Zi22[index1, , drop=FALSE]
    n <- length(index2)
    
    #���Җޓx���v�Z
    LLt21 <- LLt[index2, , , drop=FALSE]
    LLt22 <- matrix(0, nrow=n, ncol=k2)
    par <- matrix(0, nrow=n, ncol=k1)
    for(j in 1:k1){
      theta22_r <- theta22[zi22_j, , j]
      par[, j] <- (LLt21[, , j] * theta22_r) %*% vec_k2
      LLt22 <- LLt22 + LLt21[, , j] * zi21_j[, j] * theta22_r
    }
    
    #��ʊK�w�̐��ݏ�Ԑ؊����ϐ��𐶐�
    par1 <- par0 <- rep(0, n)
    for(j in 1:k1){
      par1 <- par1 + par[, j] * zi21_j[, j]
      par0 <- par0 + matrix(theta12[j, -j], nrow=n, ncol=k1-1, byrow=T) * par[, -j] * zi21_j[, j]
    }
    r <- beta[d_id[index2]]   #������
    
    #�x���k�[�C���z����؊����ϐ��𐶐�
    beta_par0 <- r * rowSums2(par0)
    beta_rate <- beta_par0 / ((1-r)*par1 + beta_par0)   #�؊����ϐ��̊����m��
    beta_rate[is.nan(beta_rate)] <- r[is.nan(beta_rate)]
    Zi1[index2] <- rbinom(n, 1, beta_rate)   
    index_z1 <- which(Zi1[index2]==1)
    
    
    ##�؊����ϐ��������Ƃ��ď�ʊK�w�𐶐�
    if(length(index_z1)==0){
      Zi21[index2, ] <- Zi21[index1, ]
      z21_vec[index2] <- z21_vec[index1]
      
    } else {
      
      r <- theta12[z21_vec[index1], ]
      Zi21[index2, ] <- Zi21[index1, ]; z21_vec[index2] <- z21_vec[index1]   #1���O�̐��ݏ�Ԃ��J��Ԃ�
      par_rate <- r*par / rowSums2(r*par)   #���ݏ�Ԃ̊����m��
      
      if(nrow(par_rate)==1){
        Zi21[index2, ] <- rmnom(length(index_z1), 1, par_rate[index_z1, ])   #�������z������ݏ�Ԃ𐶐�
        z21_vec[index2] <- as.numeric(Zi21[index2, ] %*% 1:k1)
      } else {
        Zi21[index2, ][index_z1, ] <- rmnom(length(index_z1), 1, par_rate[index_z1, ])   #�������z������ݏ�Ԃ𐶐�
        z21_vec[index2][index_z1] <- as.numeric(Zi21[index2, ][index_z1, ] %*% 1:k1)
      }
      #��ʊK�w�̃}���R�t���ڍs����X�V
      rf21 <- rf21 + t(Zi21[index1, , drop=FALSE]) %*% (Zi21[index2, , drop=FALSE] * Zi1[index2])
    }
    
    ##���ʊK�w�̐��ݏ�Ԃ𐶐�
    #��ʊK�w�̏�Ԃɉ����Ėޓx�v�Z
    par <- matrix(0, nrow=n, ncol=k2)
    for(j in 1:k1){
      index <- which(Zi21[index2, j]==1)
      index_z21 <- index[index %in% index_z1]; index_z22 <- index[!index %in% index_z1]
      if(length(index)==0) next
      LLt22 <- LLt[index2, , , drop=FALSE]
      
      #�؊����ϐ��ɉ������ޓx�v�Z
      if(length(index_z21) > 0 ){
        par[index_z21, ] <- matrix(theta21[j, ], nrow=length(index_z21), ncol=k2, byrow=T) * LLt22[index_z21, , j, drop=FALSE][, , 1]
      }
      if(length(index_z22) > 0){
        par[index_z22, ] <- theta22[z22_vec[index1][index_z22], , j] * LLt22[index_z22, , j, drop=FALSE][, , 1]
      }
    }
    
    #�������z���牺�ʊK�w�̏�Ԃ𐶐�
    par_rate <- par / rowSums2(par)   #���ݕϐ��̊����m��
    Zi22[index2, ] <- rmnom(n, 1, par_rate)   #�������z������ݏ�Ԃ𐶐� 
    z22_vec[index2] <- as.numeric(Zi22[index2, ] %*% 1:k2)
    
    #�}���R�t���ڍs����X�V
    for(j in 1:k1){
      rf22[, , j] <- rf22[, , j] + t(Zi22[index1, , drop=FALSE]) %*% (Zi22[index2, , drop=FALSE] * Zi21[index2, j] * (1-Zi1[index2]))
    }
  }
  
  ##�}���R�t���f���̃p�����[�^���T���v�����O
  #�x�[�^���z���獬�������T���v�����O
  par1 <- as.numeric(tapply(Zi1[-index_t11], d_id[-index_t11], sum))
  par2 <- w - 1 - par1
  beta <- rbeta(d, par1 + s0, par2 + v0)
  
  #�f�B���N�����z�����ʊK�w�̍��������T���v�����O
  first_vector1 <- colSums(Zi21[index_t11, ]) + alpha01
  transition_matrix1 <- rf21 + alpha01
  theta11 <- as.numeric(extraDistr::rdirichlet(1, first_vector1))
  for(j in 1:k1){
    theta12[j, -j] <- extraDistr::rdirichlet(1, transition_matrix1[j, -j])
  }
  
  #�f�B���N�����z���牺�ʊK�w�̍��������T���v�����O
  for(j in 1:k1){
    first_vector2 <- colSums(Zi22[index_t11, ] * Zi21[index_t11, j]) + colSums(Zi22 * Zi1 * Zi21[, j]) + alpha02
    transition_matrix2 <- rf22[, , j] + alpha02
    theta21[j, ] <- extraDistr::rdirichlet(1, first_vector2)
    theta22[, , j] <- extraDistr::rdirichlet(k2, transition_matrix2)
  }
  
  ##��ԂɊ�Â��g�s�b�N�𐶐�
  #�g�s�b�N�ޓx��ݒ�
  Lho_topic <- matrix(0, nrow=f, ncol=k3)
  for(j in 1:k1){
    index_z <- which(Zi21[, j]==1)
    Lho_topic[index_z, ] <- t(phi[, , j])[wd[index_z], ] * theta3[z22_vec[index_z], , j]
  }
  
  #�������z����g�s�b�N�𐶐�
  topic_rate <- Lho_topic / as.numeric(Lho_topic %*% vec_k3)   #�g�s�b�N�̊����m��
  Zi3 <- rmnom(f, 1, topic_rate)
  z3_vec <- as.numeric(Zi3 %*% 1:k3)
  
  ##�g�s�b�N���z�ƒP�ꕪ�z�𐶐�
  for(j in 1:k1){
    #�f�B���N�����z����g�s�b�N���z�𐶐�
    wsum <- t(Zi21[, j] * Zi22) %*% Zi3 + beta01
    theta3[, , j] <- extraDistr::rdirichlet(k2, wsum)
    
    #�f�B���N�����z����P�ꕪ�z�𐶐�
    vsum <- as.matrix(t(Zi21[, j] * Zi3) %*% sparse_data + beta02)
    phi[, , j] <- extraDistr::rdirichlet(k3, vsum)
  }
  
  
  ##�p�����[�^�̊i�[�ƃT���v�����O���ʂ̕\��
  #�T���v�����O���ꂽ�p�����[�^���i�[
  if(rp%%keep==0){
    #�T���v�����O���ʂ̊i�[
    mkeep <- rp/keep
    BETA[mkeep, ] <- beta
    THETA11[mkeep, ] <- theta11
    THETA12[, , mkeep] <- theta12
    THETA21[, , mkeep] <- theta21
    THETA22[, , , mkeep] <- theta22
    PHI[, , , mkeep] <- phi
  }  
  
  #�g�s�b�N�����̓o�[���C�����Ԃ𒴂�����i�[����
  if(rp%%keep==0 & rp >= burnin){
    SEG1 <- SEG1 + Zi1
    SEG21 <- SEG21 + Zi21
    SEG22 <- SEG22 + Zi22
    SEG3 <- SEG3 + Zi3
  }
  
  if(rp%%disp==0){
    #�ΐ��ޓx���v�Z
    LL <- sum(log(rowSums(Lho_topic)))
    
    #�T���v�����O���ʂ��m�F
    print(rp)
    print(c(LL, LLst))
    print(rbind(Z21_freq=colSums(Zi21), Z21_freq))
  }
}
