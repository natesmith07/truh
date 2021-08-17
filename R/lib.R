
#' Nearest neighbor computation for the TRUH statistic
#'
#' For a given \eqn{d} dimensional vector \eqn{\mathbf{y}}, this function finds the nearest neighbor of \eqn{\mathbf{y}} in
#' a \eqn{n\times d} matrix \eqn{\mathbf{U}}.
#'
#' @import Rfast
#'
#' @param y a \eqn{d} dimensional vector.
#' @param U a \eqn{n\times d} matrix where \eqn{n} represents the sample size and
#' \eqn{d} is the dimension of each sample.
#' @param n the sample size.
#' @param d dimension of each sample.
#'
#' @return
#' \enumerate{
#' \item d1 - nearest neighbor of \eqn{\mathbf{y}} in \eqn{\mathbf{U}}
#' \item d2 - nearest neighbor of d1 in \eqn{\mathbf{U}}
#' }
#'
#' @seealso \code{\link{truh}}
#'
#' @examples
#' library(truh)
#' n = 100
#' d = 3
#' set.seed(1)
#' y = rnorm(3)
#' set.seed(2)
#' U = matrix(rnorm(n*d),nrow=n,ncol=d)
#' out = nearest(y,U,n,d)
#'
#' @export

nearest<-function(y,U,n,d){
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
  return(c(d1,d2))
}

#' TRUH test statistic
#'
#' TRUH test statistic for nonparametric two sample testing under heterogeneity.
#'
#' @import Rfast
#' @import cluster
#' @import fpc
#' @import foreach
#' @import doParallel
#' @import iterators
#' @import parallel
#'
#' @param V \eqn{m\times d} matrix where \eqn{m} represents the sample size and
#' \eqn{d} is the dimension of each sample.
#' @param U a \eqn{n\times d} matrix where \eqn{n} represents the sample size and
#' \eqn{d} is the dimension of each sample with \eqn{m\ll n}.
#' @param B number of bootstrap samples.
#' @param fc fold change constant. The default value is 1. See equation (2.8) of the referenced paper for more details.
#' @param ncores the number of computing cores available. The default value is 2.
#'
#' @return
#' \enumerate{
#' \item teststat - TRUH test statistic.
#' \item k.hat - number of clusters detected in the uninfected sample.
#' \item pval - The maximum p-value across the detected clusters.
#' \item pval_all - p-value for each cluster.
#' \item dist.null_all - the approximate bootstrapped based null distribution.
#' }
#'
#' @seealso \code{\link{nearest}}
#'
#' @examples
#' library(truh)
#' n = 500
#' m = 10
#' d = 3
#' set.seed(1)
#' V = matrix(rnorm(m*d),nrow=m,ncol=d)
#' set.seed(2)
#' U = matrix(rnorm(n*d),nrow=n,ncol=d)
#' out = truh(V,U,100)
#'
#' @references
#' Banerjee, Trambak, Bhaswar B. Bhattacharya, and Gourab Mukherjee.
#' "A nearest-neighbor based nonparametric test for viral remodeling in
#' heterogeneous single-cell proteomic data."
#' The Annals of Applied Statistics 14, no. 4 (2020): 1777-1805.
#'
#' @export

truh<-function(V,U,B,fc=1,ncores=2){

  n.u<-dim(U)[1]
  n.v<-dim(V)[1]
  d.u<-dim(U)[2]
  d.v<-dim(V)[2]

  out.nearest<- matrix(0,n.v,2)
  for(i in 1:dim(V)[1]){
    out.nearest[i,]<-nearest(V[i,],U,n.u,d.u)
  }
  fac<-(n.v^{1/d.u})*(d.u>1)+1*(d.u==1)
  teststat.truh<-fac*mean(out.nearest[,1]-1*out.nearest[,2])

  #------------ cutoff computation --------------------
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  X<-U
  set.seed(1)
  k.hat<-prediction.strength(X, Gmin=2, Gmax=10,cutoff = 0.8)$optimalk
  clust.X <-clara(X,k=k.hat,metric="euclidean",samples=50,sampsize=500)$clustering
  mu<-table(clust.X)
  temp1<-diag(k.hat)
  pval<-matrix(0,k.hat,1)
  dist.null<-matrix(0,B,k.hat)

  for(k in 1:k.hat){

    s.v<-ceiling(temp1[k,k]*n.v)
    clust.idx<-which(clust.X==k)
    result<-foreach(b = 1:B,.packages="Rfast",.export="nearest")%dopar%{
      set.seed(k*b)
      idx1.v<-sample(clust.idx,s.v,replace=(s.v>mu[k]))
      VV.b<- X[idx1.v,]
      UU.b<-X[-idx1.v,]
      out<- matrix(0,dim(VV.b)[1],2)
      for(i in 1:dim(VV.b)[1]){
        out[i,]<-nearest(VV.b[i,],UU.b,dim(UU.b)[1],d.u)[1:2]
      }
      fac<-(n.v^{1/d.u})*(d.u>1)+1*(d.u==1)
      dist.null<-fac*mean(fc*out[,1]-1*out[,2])
      return(dist.null)
    }
    dist.null[,k]<-do.call(rbind,result)
    pval[k]<-sum(1*(dist.null[,k]>teststat.truh))/B
  }
  stopCluster(cl)
  registerDoSEQ()

  return(list("teststat"=teststat.truh,"pval"=max(pval),"pval_all"=pval,
              "dist.null_all"=dist.null,"k.hat"=k.hat))
}


