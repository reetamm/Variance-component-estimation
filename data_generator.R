#library(gplots)
#Settings:
#k=10,20
#N=10*k 50*k 100*k
#Imbalance= balanced, slight unbalanced, highly unbalanced



mdata<-function(k=10,n=10,sig2=0.5,sigt2=1,mu=10,mut=numeric(k)){
y<-numeric(k*n);
grp<-numeric(k*n);
out<-data.frame(y,grp);

a<-0;

for(i in 1:k){
  t=rnorm(1,mut[i],sqrt(sigt2))
for(j in 1:n){
out[j+a,1]<-mu+t+rnorm(1,0,sqrt(sig2));
out[j+a,2]<-i;
}
a<-a+n;
}
return(out);
}


mdata2<-function(k=10,n=10,b=1,sig2=0.5, m=5, sigt2=1,mu=10,mut=numeric(k)){
if(b==1){
data<-mdata(k=k,n=n,sig2=sig2,sigt2=sigt2,mu=mu,mut=mut); #changed the name of the variable from out to data
}

if(b==2){
n2<-round((k*n-2*m)/k)
out1<-mdata(k=k,n=n2,sig2=sig2,sigt2=sigt2,mu=mu,mut=mut);
out2<-mdata(k=2,n=m,sig2=sig2,sigt2=sigt2,mu=mu,mut=mut);
data<-rbind(out1,out2); #changed the name of the variable from out to data
}

if(b==3){
n2<-round((n-m)*k/2);
out1<-mdata(k=2,n=n2,sig2=sig2,sigt2=sigt2,mu=mu,mut=mut);
out2<-mdata(k=k,n=m,sig2=sig2,sigt2=sigt2,mu=mu,mut=mut);
data<-rbind(out1,out2); #changed the name of the variable from out to data
}
data$grp = as.factor(data$grp)
#plot=plotmeans(y ~ grp, data = data);
#output=list(data=data,plot=plot)

return(data)
}
#for mdata3, n_i is a vector containing the number of samples per group, k is number of groups
#by defaut, mdata3 generates balanced data with 10 groups and 10 samples per group
mdata3<-function(n_i=NULL, k=10,n=10, b=1, m=5, sig2=0.5,sigt2=1,mu=10,mut=numeric(k)){
    if(sum(is.na(n_i))==0 & sum(n_i>0)==length(n_i) & length(n_i)>0)
    {
      k<-length(n_i); #extracts number of groups
      muv<-rnorm(k,mean=0,sd=sqrt(sigt2)); #generates vector of means to be used for each group
      y<-0;
      grp<-0;
      out<-cbind(y,grp); #container of data
      for(i in 1:k){
        yd<-mu+muv[i]+rnorm(n_i[i],mean=0,sd=sqrt(sig2)); #data for group i
        out<-rbind(out,cbind(yd,numeric(n_i[i])+i)); #index i
        }
      out<-out[-1,]; #removes first row
      out<-as.data.frame(out); #turns matrix to dataframe
      out$grp = as.factor(out$grp) #coerces grp to be treated as factor
  }
else
  {
    out = mdata2(k=k,n=n,b=b,sig2=sig2,sigt2=sigt2,m=5,mu=10,mut=numeric(k))
  }
#plotmeans(y ~ grp, data = out); #plot
return(out);
}

mdata3(n_i=c(NA))
