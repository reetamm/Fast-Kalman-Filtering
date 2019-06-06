#install.packages("expm")
library(Matrix)
library(expm)
library(MASS)
library(matrixcalc)
library(ks)
library(ggplot2)
library(tidyr)

rmse = function(x)
{
  n = length(x)
  out = sqrt((1/n)*sum(x^2))
  return(out)
}

T = 500 #total sample size
d = 100  #dimension of state vector
b = 1   #dimenstion of y
A = 0.95 * diag(d)  #State matrix dxd
W = as.matrix(0.5) #Variance of observation error bxb
V = 0.1 * diag(d) #Variance of state error dxd
B = matrix(rep(b,d),nrow = 1) #observation gain matrix bxd


### Data generation

#vecV0 = solve(diag(d**2) - A %x% A) %*% vec(V)
#V0 = invvec(vecV0)

V0 = diag(d) - diag(d)

for(i in 1:1000)
{
  status = paste('iteration',i,', V0[1,1] =',V0[1,1])  
  print(status)
  tol = V0[1,1]
  V0 = V0 + (A %^% i) %*% V %*% t(A %^% i)
  if(abs(tol-V0[1,1])<.000001)
    break
}
V0[1,1]

x0 = mvrnorm(1,mu = rep(0,d),Sigma = V0)
x = matrix(0, nrow = d, ncol = T)
y = numeric(T)
x[,1] = A %*% x0 + mvrnorm(1,rep(0,d),V)
y[1] = B %*% x0 + mvrnorm(1,rep(0,b),W)
for(i in 2:T)
{
  print(i)
  x[,i] = A %*% x[,i-1] + mvrnorm(1,rep(0,d),V)
  y[i] = B %*% x[,i] + mvrnorm(1,rep(0,b),W)
}
#########################################################
########################################################
#x50 = x
#y50=y
#x100=x
#y100=y
#x250=x
#y250=y
#x500=x
#y500=y
#x1000 = x
#y1000 = y
####################Fast Kalman Filtering ###############
x=x1000
y=y1000
mu.fkf = matrix(nrow = d, ncol = T)
mu.tminus1 = rep(0,d)


start.fkf = proc.time()
#P.t = vector('list')
C0 = V %*% solve(diag(d) - A %*% t(A))  #base for the low rank approx of C_t, dxd
C0invA = solve(C0) %*% A   #inverse(C0) x A, dxd

L.t = C0 %*% t(B)  #dxb
SIGinv = W + B %*% L.t #bxb
SIG = solve(SIGinv)
mu.fkf[,1] = L.t %*% SIGinv %*% y[1]  #conditional mean of x_t given y(1:t), dx1
#P.t[[1]] = C0 - L.t %*% SIG %*% t(L.t)  # dxd

for(t in 2:T)
{
  #print(t) 
  mu.tminus1 = mu.fkf[,t-1]
  
  AL = A %*% L.t  # dxb
  Phi = C0invA %*% L.t  # dxb intermediate matrix
  Dinv = SIGinv - t(AL) %*% Phi # bxb intermediate matrix
  O = cbind(Phi,t(B)) # dx2b intermediate matrix
  Qinv = bdiag(Dinv,W)  # 2bx2b
  decomp = svd(C0 %*% O %*% solve(sqrtm(Qinv + t(O) %*% C0 %*% O))) #svd dx2b
  Lh = decomp$u # dxb
  #SIGh = diag(decomp$d) %^% 2 #2bx2b
  SIGh = decomp$d^2
  L.t = as.matrix(Lh[,1])
  
  P = C0 - AL %*% SIG %*% t(AL)
  mu.fkf[,t] = A %*% mu.tminus1 + P %*% t(B) %*% solve(W + B %*% P %*% t(B)) %*% (y[t] - B %*% A %*% mu.tminus1)

  SIG = SIGh[1]
  SIGinv = solve(SIG)
  #P.t[[t]] = C0 - L.t %*% SIG %*% t(L.t)
  
}
end.fkf = proc.time()
time.fkf = end.fkf - start.fkf
time.fkf
##############
#Kalman filter
mu.kf = matrix(nrow = d, ncol = T)

start.kf = proc.time()
C0 = V %*% solve(diag(d) - A %*% t(A))  #base for the low rank approx of C_t, dxd
for(t in 1:T)
{
  print(t)
  if(t==1)
  {
C.pr = C0
mu.tminus1 = rep(0,d)
  }
  else
  {
    C.pr = C.t
    mu.tminus1 = mu.kf[,t-1]
  }
  P = A %*% C.pr %*% t(A) + V
  PtB = P %*% t(B)
  K.t = PtB %*% solve(B %*% PtB + W) 
  C.t = P -  K.t %*% t(PtB)
  mu.kf[,t] = A %*% mu.tminus1 + K.t %*% (y[t] - B %*% A %*% mu.tminus1)
}
end.kf = proc.time()
time.kf = end.kf - start.kf

time.fkf ; time.kf

#########################
#performance

#mu50.kf = mu.kf
#mu50.fkf = mu.fkf
yhat50.fkf = apply(mu50.fkf,2,sum)

mu100.kf = mu.kf
mu100.fkf = mu.fkf
yhat100.fkf = apply(mu100.fkf,2,sum)
yhat100.kf = apply(mu100.kf,2,sum)


mu250.kf = mu.kf
mu250.fkf = mu.fkf
yhat250.fkf = apply(mu250.fkf,2,sum)
yhat250.kf = apply(mu250.kf,2,sum)

#mu1000.kf = mu.kf
#mu1000.fkf = mu.fkf
yhat1000.fkf = apply(mu1000.fkf,2,sum)

rmse1000 = apply(mu1000.fkf-x1000,2,rmse)
rmse250 = apply(mu250.fkf-x250,2,rmse)
rmse100 = apply(mu100.fkf-x100,2,rmse)
rmse50 = apply(mu50.fkf-x50,2,rmse)
rmse250.kf = apply(mu250.kf-x250,2,rmse)
rmse100.kf = apply(mu100.kf-x100,2,rmse)

rmsefkf = data.frame(time=4:T,d50 = rmse50[4:T],d100 = rmse100[4:T],d250 = rmse250[4:T],d1000 = rmse1000[4:T])
rmselong = gather(rmsefkf,d,rmse,d50:d1000,factor_key=T)
ggplot(rmselong,aes(x=time,y=rmse,col=d)) + geom_line() + scale_color_manual(values = c('red','blue','green','black'))
ggplot(rmselong,aes(x=d,y=rmse,fill=d)) + geom_boxplot()

apply(rmsefkf, 2, sum)

rmsecomp = data.frame(time=4:T,kf=rmse250.kf[4:T],fkf=rmse250[4:T])
rmsecomplong = gather(rmsecomp,method,rmse,kf:fkf,factor_key = T)
ggplot(rmsecomplong,aes(x=time,y=rmse,col=method)) + geom_line() + scale_color_manual(values = c('red','black')) + labs(title = 'RMSE for KF and FKF when d=250')
sum(rmse50)
