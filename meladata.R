preprocessed<-function(data){
  len=length(data[,1])
  c1=matrix(0,nrow=len,ncol=1)
  c2=matrix(0,nrow=len,ncol=1)
  c3=matrix(0,nrow=len,ncol=1)
  c4=matrix(0,nrow=len,ncol=1)
  for(i in 1:len){
    if(data$x[i]==1){
      c1[i]=1
    }
    if(data$x[i]==2){
      c2[i]=1
    }
    if(data$x[i]==3){
      c3[i]=1
    }
    if(data$x[i]==4){
      c4[i]=1
    }
  }
  data=cbind(data,c1,c2,c3,c4)
  FAVOR=matrix(0,nrow=len,ncol=1)
  for(i in 1:len){
    
    if(i>=(len-10)){
      FAVOR[i]=1
    }
    if(i<(len-10)&i>=10){
      z=runif(1,0,1)
      if(z>=i/len){
        FAVOR[i]=0
      }
      else{
        FAVOR[i]=1
      }
    }
  }
  data=cbind(data,FAVOR)
  CURED=matrix(0,nrow=len,ncol=1)
  for(i in 1:len){
    if(data$d[i]==1){
      CURED[i]=0
    }
    else{
      s=runif(1,0,1)
      if(s>=i/length(data[,1])){
        CURED[i]=0
      }
      else{
        CURED[i]=1
      }
    }
  }
  data=cbind(data,CURED)
  return(data)
}
processed<-function(data){
  len=length(data[,1])
  a=subset(data,data$d==1)
  b=a[!duplicated(a$t),]
  dup=length(a[,1])
  nodup=length(b[,1])
  interval=matrix(0,nrow=len,ncol=1)
  for(i in 1:len){#compute interval
    time=data$t[i]
    for(j in 1:(nodup-1)){
      if(time>=b$t[j]&time<b$t[j+1]){
        interval[i]=j  
        break
      }}
    if(time>=b$t[nodup]){
      interval[i]=j+1
    }
  }
  data=cbind(data,interval)
  count=matrix(0,nrow=nodup,ncol=1)
  for(i in 1:nodup){  
    time=b$t[i]
    count[i]=length(which(data$t*data$d==time))
  }
  return(list(data,count))
}
m1234=function(data,i,delta,pi1,pi2,lambda2,alpha,beta2,s_t,z)
{
  s0=s_t[data$interval[i]]
  s1=s0^exp(beta2%*%z[i,2:4])
  s2=s0^exp(beta2%*%z[i,2:4]+lambda2)
  p=(s2*pi2+1-pi2)*(1-alpha)+(s1*pi1+1-pi1)*alpha
  m1=delta*pi2*s2*(1-alpha)/p+
    (1-delta)*s2*exp(lambda2)*(1-alpha)*pi2/(s2*exp(lambda2)*(1-alpha)*pi2+s1*alpha*pi1)
  m2=delta*(1-pi2)*(1-alpha)/p
  m3=delta*pi1*s1*alpha/p+(1-delta)*alpha*pi1*s1/(alpha*pi1*s1+s2*exp(lambda2)*(1-alpha)*pi2)
  m4=delta*(1-pi1)*alpha/p
  
  return(c(m1,m2,m3,m4))
}

st<-function(data,beta2,lambda2,z,count,m)
{
  #Survival function at time t
  a=numeric()
  s=numeric()
  n=dim(data)[1]
  for(i in 1:max(data$interval))
  {
    first=min(which(data$interval==i))
    zz=z[first:n,]
    mm=m[first:n,]
    a[i]=-count[i]/sum((mm[,1]*exp(lambda2)+mm[,3])*exp(zz[,2:4]%*%beta2))
    s[i]=exp(sum(a))
  }
  
  return(s)
}

ht<-function(data,beta2,lambda2,z,count,m)
{
  #Hazard function at time t
  a=numeric()
  s=numeric()
  n=dim(data)[1]
  for(i in 1:max(data$interval))
  {
    first=min(which(data$interval==i))
    zz=z[first:n,]
    mm=m[first:n,]
    a[i]=count[i]/sum((mm[,1]*exp(lambda2)+mm[,3])*exp(zz[,2:4]%*%beta2))
  }
  return(a)
}

# M-step
L1<-function(data,par,m,z,delta)
{
  pi1=exp(z%*%as.vector(c(par[1],par[2],par[3],par[4])))    
  pi2=pi1*exp(par[5])    
  pi1=1/(1+pi1)
  pi2=1/(1+pi2)
  
  L1=sum(m[,1]*log(pi2)+m[,2]*log(1-pi2)+m[,3]*log(pi1)+m[,4]*log(1-pi1)) 
  
  return(-L1)
}

L2<-function(data,m,d,pa,count,z)
{ 
  L2=0
  n=dim(data)[1]
  for(i in 1:length(count))
  {
    k=min(which(data$interval==i))
    z1=z[k:n,2]
    z2=z[k:n,3]
    z3=z[k:n,4]
    m1=m[k:n,1] 
    m3=m[k:n,3]
    summ=sum(m1*exp(pa[1]*z1+pa[2]*z2+pa[3]*z3+pa[4])+m3*exp(pa[1]*z1+pa[2]*z2+pa[3]*z3))
    L2=L2+log(exp(pa[1]*sum(z1[1:count[i]])+pa[2]*sum(z2[1:count[i]])+pa[3]*sum(z3[1:count[i]])+
                    sum(d[k:(k+count[i]-1)])*pa[4])/(summ)^(count[i]))
    
  }
  return(-L2)
}

L3<-function(data,d,pa,z)
{
  alpha=exp(z%*%as.vector(c(pa[1],pa[2],pa[3],pa[4])))
  alpha=1/(1+alpha)
  L3=sum(d*log(1-alpha)+(1-d)*log(alpha))
  return(-L3)
}
library(survival)
library(smcure)
library(dplyr)
data1=read.csv(file = "melanomadata.csv")
data1=data1%>%arrange(t)
data2=preprocessed(data1)
data3=processed(data2)
zero_1=data3[[1]]
count=data3[[2]]
len=length(zero_1[,1])

logi1=glm(CURED~c1+c2+c3+c4+FAVOR,family = binomial(link = "logit"),data=zero_1)
beta1<-logi1$coefficients[1:4]
lambda1<-logi1$coefficients[6]

sur<-Surv(zero_1$t,zero_1$d)
logi2=coxph(sur~c1+c2+c3+c4+FAVOR,data=zero_1)
beta2<-logi2$coef[1:3]
lambda2<-logi2$coef[5]

logi3=glm(FAVOR~c1+c2+c3+c4,family = binomial(link = "logit"),data=zero_1)
w<-logi3$coef[1:4]

p<-matrix(0,dim(zero_1)[1],ncol=1)
m<-matrix(0,dim(zero_1)[1],ncol=4)
for(i in 1:dim(zero_1)[1])
{
  if(zero_1$FAVOR[i]==1 & zero_1$CURED[i]==0){m[i,1]=1}
  if(zero_1$FAVOR[i]==1 & zero_1$CURED[i]==1){m[i,2]=1}
  if(zero_1$FAVOR[i]==0 & zero_1$CURED[i]==0){m[i,3]=1}
  if(zero_1$FAVOR[i]==0 & zero_1$CURED[i]==1){m[i,4]=1}
}

d<-matrix(0,dim(zero_1)[1],ncol=1)

error=10
times=0
delta=1-zero_1$d
z=as.matrix(cbind(1,zero_1[,"c1"],zero_1[,"c2"],zero_1[,"c3"]))


while(error>10^(-6))
{ 
  pi1=exp(z%*%beta1)
  pi2=1/(1+pi1*exp(lambda1))
  pi1=1/(1+pi1)
  alpha=1/(1+exp(z%*%w))   #1 unfavor 
  s_t=st(zero_1,beta2,lambda2,z,count,m)
  
  #print(m)
  for(i in 1:dim(zero_1)[1])
  { 
    m[i,]=m1234(zero_1,i,delta[i],pi1[i],pi2[i],lambda2,alpha[i],beta2,s_t,z)
    d[i]=m[i,1]+m[i,2]
  }
  
  #options(max.print=1000000)  
  #print(m)
  #print(zero_1)
  result1=optim(par=c(beta1,lambda1),fn=L1,data=zero_1,m=m,z=z,delta=delta)
  result2=optim(par=c(beta2,lambda2),fn=L2,data=zero_1,m=m,d=d,count=count,z=z)
  result3=optim(par=w,fn=L3,data=zero_1,d=d,z=z)
  
  beta1_new=result1$par[1:4]
  lambda1_new=result1$par[5]
  
  beta2_new=result2$par[1:3]
  lambda2_new=result2$par[4]
  w_new=result3$par
  
  
  error2=max(abs(lambda1_new/lambda1-1))
  error1=max(abs(beta1_new/beta1-1))
  error3=max(abs(beta2_new/beta2-1))
  error4=max(abs(lambda2_new/lambda2-1))
  error5=max(abs(w_new/w-1))
  error=max(c(error1,error2,error3,error4,error5))
  
  beta1=beta1_new
  lambda1=lambda1_new
  beta2=beta2_new
  lambda2=lambda2_new
  w=w_new
  
  
  times=times+1
  
  print(times)
  print(error)
}

pi1=exp(z%*%beta1)
pi2=1/(1+exp(z%*%beta1+lambda1))
pi1=1/(1+pi1)
alpha=1/(1+exp(z%*%w)) 
s_t=st(zero_1,beta2,lambda2,z,count,m) 
h_t=ht(zero_1,beta2,lambda2,z,count,m) 

for(i in 1:dim(zero_1)[1])
{ 
  s0=s_t[zero_1$interval[i]]
  h0=h_t[zero_1$interval[i]]
  s1=s0^exp(beta2%*%z[i,2:3])
  s2=s0^exp(beta2%*%z[i,2:3]+lambda2)
  
  p[i]=log((1-alpha[i])*(pi2[i]*h0*s2*exp(beta2%*%z[i,2:3]+lambda2))^(1-delta[i])*(s2*pi2[i]+1-pi2[i])^delta[i]+
             alpha[i]*(pi1[i]*h0*s1*exp(beta2%*%z[i,2:3]))^(1-delta[i])*(s1*pi1[i]+1-pi1[i])^delta[i])
} 
print(sum(p))
