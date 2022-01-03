KM<-function(data){
  #making interval for censoring
  len=length(data[,1])
  a=subset(data,data$d==0)
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
  datakm=cbind(data,interval)
  
  d=matrix(0,nrow=max(datakm$interval),ncol=1) #uncensored+censored
  r=matrix(0,nrow=max(datakm$interval),ncol=1) #risk set
  G=matrix(0,nrow=max(datakm$interval),ncol=1) #Kaplan-Meier estimate
  gp=matrix(0,nrow=max(datakm$interval),ncol=1)
  for(i in 1:max(datakm$interval)){ 
    d[i]=length(which(datakm$interval==i))
  }
  r[1]=dim(datakm)[1]-length(which(datakm$interval==0))
  for(i in 2:max(datakm$interval)){ 
    r[i]=dim(datakm)[1]-length(which(datakm$interval==0))-sum(d[1:i-1,])
  }
  for(i in 1:max(datakm$interval)){
    G[i]=(r[i]-1)/r[i] #s[i]==1
    gp[i]=prod(G[1:i,])
  }
  gp=1-gp
  return(list(datakm,gp))
}
modeling1<-function(data,tRbeta1,tRlambda1,tRbeta2,tRlambda2,tRw){
  data3=processed(data%>%arrange(t))
  zero_1=data3[[1]]
  count=data3[[2]]
  len=length(zero_1[,1])
  
  beta1=tRbeta1
  lambda1=tRlambda1
  beta2=tRbeta2
  lambda2=tRlambda2
  w=tRw
  
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
    
    
    #times=times+1
    
    #print(times)
    #print(error)
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
    s1=s0^exp(beta2%*%z[i,2:4])
    s2=s0^exp(beta2%*%z[i,2:4]+lambda2)
    
    p[i]=log((1-alpha[i])*(pi2[i]*h0*s2*exp(beta2%*%z[i,2:4]+lambda2))^(1-delta[i])*(s2*pi2[i]+1-pi2[i])^delta[i]+
               alpha[i]*(pi1[i]*h0*s1*exp(beta2%*%z[i,2:4]))^(1-delta[i])*(s1*pi1[i]+1-pi1[i])^delta[i])
  } 
  print(sum(p))
  return(list(beta1,lambda1,beta2,lambda2,w,zero_1,count,s_t,pi1,pi2,alpha,z))
}
bootstrap <- function(zero_1,count,s_t,pi1,pi2,alpha,beta2,z,lambda2,datakm,km){
  data_1=zero_1
  #Step 1 and 2
  for(i in 1:dim(zero_1)[1]){
    data_1$FAVOR[i]=rbinom(1,1,1-alpha[i]) 
    if(data_1$FAVOR[i]==0){
      data_1$CURED[i]=rbinom(1,1,1-pi1[i])
    }
    else{
      data_1$CURED[i]=rbinom(1,1,1-pi2[i])
    }
  }
  #Step 3 
  for(i in which(data_1$CURED==0)){
    if(data_1$FAVOR[i]==0){
      F=1-s_t^as.numeric(exp(beta2 %*% z[i, 2:4]))
    }
    else{
      F=1-s_t^as.numeric(exp(beta2 %*% z[i, 2:4]+lambda2))
    }
    rv=runif(1,0,1)
    if(length(which(F<rv))<length(count)){
      TF=zero_1$t[min(which(zero_1$interval==((length(which(F<rv)))+1)))]
    }
    else{
      TF=zero_1$t[min(which(zero_1$interval==length(count)))]
    }
    
    #Step 4
    SG=matrix(0,nrow=max(datakm$interval),ncol=1)
    if(datakm$d[i]==0){
      TC=datakm$t[i]
    }
    
    else{
      if(datakm$interval[i]<max(datakm$interval)){
        if(datakm$interval[i]==0){
          for(k in (datakm$interval[i]+1):max(datakm$interval)){
            SG[k]=km[k]
          }
        }
        else{
          for(k in (datakm$interval[i]+1):max(datakm$interval)){
            SG[k]=(km[k]-km[datakm$interval[i]])/(1-km[datakm$interval[i]])
          }
        }
        
        rn=runif(1,0.1)
        if(length(which(SG<rn))<max(datakm$interval)){
          TC=datakm$t[min(which(datakm$interval==((length(which(SG<rn)))+1)))]
        }
        else{
          TC=datakm$t[min(which(datakm$interval==max(datakm$interval)))]
        }
      }
      
      else{
        TC=datakm$t[min(which(datakm$interval==max(datakm$interval)))]
      }
      
    }
    data_1$t[i]=min(TC,TF)
    data_1$d[i]=as.numeric(TC>TF)
    #if TC=TF set delta=0.
  }
  for(i in which(data_1$CURED==1)){
    data_1$d[i]=0
    ru=runif(1,0,1)
    if(length(which(km<ru))<max(datakm$interval)){
      TC=datakm$t[min(which(datakm$interval==((length(which(km<ru)))+1)))]
    }
    else{
      TC=datakm$t[min(which(datakm$interval==max(datakm$interval)))]
    }
  }
  data_1=subset(data_1,select=-interval)
  return(data_1)
  
}

library(survival)
library(smcure)
library(dplyr)

tRzero_1=lista[[6]]
tRcount=lista[[7]]
tRs_t=lista[[8]]
tRpi1=lista[[9]]
tRpi2=lista[[10]]
tRalpha=lista[[11]]
tRz=lista[[12]]
tRbeta2=lista[[3]]
tRlambda2=lista[[4]]
tRbeta1=lista[[1]]
tRlambda1=lista[[2]]
tRw=lista[[5]]
data_km=subset(tRzero_1,select=-interval)
listb=KM(data_km)
datakm=listb[[1]]
km=listb[[2]]

Beta11=numeric()
Beta12=numeric()
Beta13=numeric()
Lambda1=numeric()
Beta21=numeric()
Beta22=numeric()
Lambda2=numeric()
W1=numeric()
W2=numeric()
W3=numeric()
timess=0
for(i in 1:5){
  datat=bootstrap(tRzero_1,tRcount,tRs_t,tRpi1,tRpi2,tRalpha,tRbeta2,tRz,tRlambda2,datakm,km)
  listbs=tryCatch(modeling1(datat,tRbeta1,tRlambda1,tRbeta2,tRlambda2,tRw),error=function(e) NULL,warning=function(e) NULL)
  if(class(listbs)=="list"){
    Beta11[i]=listbs[[1]][1]
    Beta12[i]=listbs[[1]][2]
    Beta13[i]=listbs[[1]][3]
    Lambda1[i]=listbs[[2]]
    Beta21[i]=listbs[[3]][1]
    Beta22[i]=listbs[[3]][2]
    Lambda2[i]=listbs[[4]]
    W1[i]=listbs[[5]][1]
    W2[i]=listbs[[5]][2]
    W3[i]=listbs[[5]][3]
  }
  timess=timess+1
  print(timess)
  
}
print(Beta11)
print(Beta12)
print(Beta13)
print(Lambda1)
print(Beta21)
print(Beta22)
print(Lambda2)
print(W1)
print(W2)
print(W3)

#Beta11=na.omit(Beta11)


#quantile(Beta11,c(0.025,0.975))

