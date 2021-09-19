library(readr)
dvis <- read_csv("D:/ѧϰ/MSc/Practicals/CompStats-Week1-TT21/dvis.csv")

glm <- glm(docvis~1+age+hhninc+female+hhkids+educyrs+addins,
           data=dvis, family=poisson())
beta.hat <- glm$coefficients
summary(glm)

######################################################
# (Intercept)  age      hhninc   female   hhkids     #
# 0.182105     0.002273 0.015997 0.048738 0.053270   #
# educyrs      addins                                #
# 0.010740     0.126796                              # 
######################################################
  
library(sandwich)
rcov1 <- vcovHC(glm,type='HC0')
se1 <- sqrt(diag(rcov1))
se1

###############################################################
# (Intercept)     age      hhninc      female      hhkids     #
# 0.254294467 0.003312846 0.020792797 0.068652805 0.074212943 #
# educyrs      addins                                         #
# 0.014424773 0.194750943                                     # 
###############################################################

B <- 1e5
n <- nrow(dvis)
beta.mat <- matrix(data=NA,nrow=B,ncol=7)
for (i in 1:B){
  rows <- sample(1:n, replace=T)
  glm.boot <- glm(docvis~1+age+hhninc+female+hhkids+educyrs+addins,
              data=dvis[rows,], family=poisson())
  beta.boot <- glm.boot$coefficients
  beta.mat[i,] <- beta.boot
  if (i%%4000==0) print('4000 boots finished')
}

se.boot <- sqrt(diag(var(beta.mat)))
###############################################################
# (Intercept)         age      hhninc      female      hhkids #
# 0.255553659 0.003329749 0.020978023 0.068824701 0.074631902 #
# educyrs      addins                                         #
# 0.014497522 0.205093932                                     # 
###############################################################

betas <- beta.hat[c('hhkids','educyrs')]
W.glm <- t(betas)%*%solve(vcov(glm)[5:6,5:6])%*%betas #8.71, p=0.09216
W.sandwich <-  t(betas)%*%solve(rcov1[5:6,5:6])%*%betas #4.86, p=0.08803683
W.boot <- t(betas)%*%solve(var(beta.mat)[5:6,5:6])%*%betas #4.79, p=0.09117268

W.vec <- vector(length=B)
for (i in 1:B){
  rows <- sample(1:n, replace=T)
  glm.boot <- glm(docvis~1+age+hhninc+female+hhkids+educyrs+addins,
                  data=dvis[rows,], family=poisson())
  beta.boot <- glm.boot$coefficients
  beta.s <- beta.boot[c('hhkids','educyrs')]
  rcov.boot <- vcovHC(glm.boot,type='HC0')
  W.vec[i] = t(beta.s-beta.hat[c('hhkids','educyrs')])%*%
    solve(rcov.boot[5:6,5:6])%*%(beta.s-beta.hat[c('hhkids','educyrs')])
  
  if (i%%4000==0) print('4000 boots finished')
}

p.value <- sum(W.vec>4.86)/B


W.vec.para <- vector(length=B)
glm0 <- glm(docvis~1+age+hhninc+female+addins,data=dvis, family=poisson())
mu.hat0 <- glm0$fitted.values
for (i in 1:B){
  y.boot <- rpois(1204, mu.hat0)
  glm.boot.para <- glm(y.boot~1+age+hhninc+female+hhkids+educyrs+addins,
                  data=dvis, family=poisson())
  beta.boot.para <- glm.boot.para$coefficients
  beta.s <- beta.boot.para[c('hhkids','educyrs')]
  rcov.boot <- vcovHC(glm.boot.para,type='HC0')
  W.vec.para[i] = t(beta.s)%*%solve(rcov.boot[5:6,5:6])%*%beta.s
  
  if (i%%4000==0) print('4000 boots finished')
}

p.value.para <- sum(W.vec.para>4.86)/B
