library(Rfast)
library(ggplot2)
library(fpc)
library(cluster)

nearest<-function(y,U,n=100,d=10){
  b1=matrix(y,ncol=d,nrow=n,byrow=T)
  temp=rowsums(abs(b1-U))
  ind1=which.min(temp);
  y.u=U[ind1,];
  d1=temp[ind1];
  
  b2=matrix(y.u,ncol=d,nrow=n-1,byrow=T)
  U1=as.matrix(U[-c(ind1),]);
  temp=rowsums(abs(b2-U1))
  ind2=which.min(temp);
  u.u=U1[ind2,];
  d2=temp[ind2];
  return(c(d1,d2,d1-d2))
}

nearest.nulldist<-function(y,U,n=100,d=10){
  b1=matrix(y,ncol=d,nrow=n,byrow=T)
  temp=rowsums(abs(b1-U))
  ind1=which.min(temp);
  y.u=U[ind1,];
  d1=temp[ind1];
  
  b2=matrix(y.u,ncol=d,nrow=n-1,byrow=T)
  U1=as.matrix(U[-c(ind1),]);
  temp=rowsums(abs(b2-U1))
  ind2=which.min(temp[temp!=0]);
  d2=temp[ind2];
  return(c(d1,d2,d1-d2))
}

truh<-function(V,U,B,fc=1,a=0.05){
  
  n.u<-dim(U)[1]
  n.v<-dim(V)[1]
  d.u<-dim(U)[2]
  d.v<-dim(V)[2]
  
  out.nearest<- matrix(0,n.v,3)
  for(i in 1:dim(V)[1]){
    out.nearest[i,]<-nearest(V[i,],U,n.u,d.u)
  }
  fac<-(n.v^{1/d.u})*(d.u>1)+1*(d.u==1)
  teststat.truh<-fac*mean(out.nearest[,1]-1*out.nearest[,2])
  
  #------------ cutoff computation --------------------
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  X<-U
  set.seed(1)
  k.hat<-prediction.strength(X, Gmin=2, Gmax=10,cutoff = 0.8)$optimalk
  clust.X <-clara(X,k=k.hat,metric="euclidean",samples=50,sampsize=500)$clustering
  mu<-table(clust.X)
  temp1<-diag(k.hat)
  w.v.b<-temp1[rep(1:k.hat,each=B),]
  s.v<- ceiling(w.v.b*n.v)
  Bp<-dim(s.v)[1]
  if(k.hat==1){
    Bp = B
  }
  
  result<-foreach(b = 1:Bp,.packages="Rfast",.export="nearest")%dopar%{
    
    list.u.b<-list()
    list.v.b<-list()
    if(k.hat>1){
      for(k in 1:k.hat){
        clust.idx<-which(clust.X==k)
        
        set.seed(k*b)
        idx1.v<-sample(clust.idx,s.v[b,k],replace=(s.v[b,k]>mu[k]))
        idx1.u<-clust.idx[!(clust.idx %in% idx1.v)]
        list.v.b[[k]]<- X[idx1.v,]
        list.u.b[[k]]<-X[idx1.u,]
      }
    } else{
      clust.idx<-which(clust.X==1)
      
      set.seed(b)
      idx1.v<-sample(clust.idx,s.v[b],replace=(s.v[b]>mu[1]))
      idx1.u<-clust.idx[!(clust.idx %in% idx1.v)]
      list.v.b[[1]]<- X[idx1.v,]
      list.u.b[[1]]<-X[idx1.u,]
      
    }
    if(d.u>1){
    VV.b<- do.call(rbind,list.v.b)
    UU.b<- do.call(rbind,list.u.b)
    } else {
      VV.b<-as.matrix(unlist(list.v.b))
      UU.b<-as.matrix(unlist(list.u.b))
    }
    out<- matrix(0,dim(VV.b)[1],2)
    for(i in 1:dim(VV.b)[1]){
      out[i,]<-nearest(VV.b[i,],UU.b,dim(UU.b)[1],d.u)[1:2]
    }
    
    fac<-(n.v^{1/d.u})*(d.u>1)+1*(d.u==1)
    dist.null<-fac*mean(fc*out[,1]-1*out[,2])
  }
  dist.null<-do.call(rbind,result)
  #*******************************************************************
  c1<-quantile(dist.null,a/2)
  c2<-quantile(dist.null,1-a/2)
  pval<-sum(1*(dist.null>teststat.truh))/Bp
  stopCluster(cl)
  registerDoSEQ()
  return(list("teststat"=teststat.truh,"c1"=c1,"c2"=c2,"pval"=pval,"null.dist"=dist.null))
  
}
