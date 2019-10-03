# Variance of rows for arrays with NA values
f.rvar.na <- function(x) {
  ave <- as.vector(f.rmean.na(x))
  pos <- which(is.na(x))
  largo <- f.rlength.na(x)
  x[pos] <- rep(ave,ncol(x))[pos]
  (x-ave)^2 %*% rep(1,ncol(x))/(largo-1)
}

# Mean of rows for arrays with NA values
f.rmean.na <- function(x) {
  largo <- f.rlength.na(x)
  x[is.na(x)] <- 0
  res <- x %*% rep(1,ncol(x)) / largo
  names(res) <- rownames(x)
  res
}

# Length of rows for arrays with NA values
f.rlength.na <- function(x) {
  r <- x/x
  r[x==0] <- 1
  r[!is.finite(r)] <- 0
  r %*% rep(1,ncol(r))
}

# Function to extract correlation indeces
extract<-function(mat,x,y){
  mat[x+nrow(mat)*(y-1)]
}

# Function to process combos
funcombos<-function(inp,tmat,groups,r){

  i<-inp[1]
  j<-inp[2]
  submat1<-tmat[,groups[[i]]]
  submat2<-tmat[,groups[[j]]]
  cormat<-cor(submat1,submat2)

  # Fix problem with symmetrical quadrants
  if(i==j){
    cormat[lower.tri(cormat,diag=TRUE)]<-0
  }

  # Extract significant correlations
  hits<-which(abs(cormat)>=r,arr.ind=TRUE)
  rowhits<-rownames(cormat)[hits[,1]]
  colhits<-colnames(cormat)[hits[,2]]

  # Extract correlation indeces
  corhits<-cormat[hits[,1]+nrow(cormat)*(hits[,2]-1)]

  # Results
  results<-as.data.frame(cbind(rowhits,colhits,corhits))

  return(results)
}


# # Quadrant Fast correlation
# qfcor<-function(inmat,ncores=1){
#   tmat<-t(inmat)
#   nfeatures<-ncol(tmat)
#   features<-colnames(tmat)
#
#   # Split into m correlation matrices of size n
#   n<-1000
#   groups<-split(features,ceiling(seq_along(features)/n))
#   m<-length(groups)
#
#   # Create quadrants combinations
#   combos<-list()
#   ij<-1
#   for (i in 1:length(groups)){
#     for(j in i:length(groups)){
#       combos[[ij]]<-c(i,j)
#       ij<-ij+1
#     }
#   }
#
#   # Calculate m x m correlation matrices in parallel
#   cl<-makeCluster(ncores)
#   siglist<-clusterApply(cl,combos,funcombos,tmat=tmat,groups=groups,r=r)
#   stopCluster(cl)
#
#   # Merge results into a final table
#   sigedges<-matrix(nrow=0,ncol=3)
#   for(i in 1:length(siglist)){
#     sigedges<-rbind(sigedges,siglist[[i]])
#   }
#
#   return(sigedges)
# }


# Fast correlation
fcor<-function(inmat,centroids,r){
  tmat<-t(inmat)
  nfeatures<-ncol(tmat)
  features<-colnames(tmat)
  targets<-setdiff(features,centroids)

  # Calculate centroid x target correlations
  cenmat<-tmat[,centroids]
  tarmat<-tmat[,targets]
  cormat<-cor(cenmat,tarmat)

  # Extract significant correlations
  hits<-which(abs(cormat)>=r,arr.ind=TRUE)
  rowhits<-rownames(cormat)[hits[,1]]
  colhits<-colnames(cormat)[hits[,2]]

  # Extract correlation indeces
  corhits<-cormat[hits[,1]+nrow(cormat)*(hits[,2]-1)]

  # Results
  sigedges<-as.data.frame(cbind(rowhits,colhits,corhits))

  return(sigedges)
}

# Function to bootstrap
bootfun<-function(inmat,centroids,r,seed=NULL){
  set.seed(seed)
  bootmat<-inmat[,sample(colnames(inmat),replace=TRUE)]
  # Calculate correlations in the bootstrapped matrix
  options(warn=-1)
  bootsigedges<-fcor(bootmat,centroids,r)
  options(warn=0)
  return(bootsigedges)
}

