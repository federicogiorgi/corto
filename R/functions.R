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

# Max of absolute value
maxabs<-function(x){
  max(abs(x))
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




#' A fast correlation function
#' @param inmat An input matrix with features as rows and samples as columns
#' @param centroids A character vector indicating the centroids
#' @param r A numeric correlation threshold
#' @return A matrix describing which edges were significant in the input matrix
#' matrix according to the r correlation threshold provided
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

#' Function to bootstrap matrix
#' @param inmat An input matrix with features as rows and samples as columns
#' @param centroids A character vector indicating the centroids
#' @param r A numeric correlation threshold
#' @param seed An integer to set the random seed
#' @return A matrix describing which edges were significant in the bootstrapped
#' matrix according to the r correlation threshold provided
bootmat<-function(inmat,centroids,r,seed=NULL){
  set.seed(seed)
  bootmat<-inmat[,sample(colnames(inmat),replace=TRUE)]
  # Calculate correlations in the bootstrapped matrix
  options(warn=-1)
  bootsigedges<-fcor(bootmat,centroids,r)
  options(warn=0)
  return(bootsigedges)
}

#' Function to process bootstraps into edges
#' @param seed An integer to set the random seed
#' @param inmat An input matrix with features as rows and samples as columns
#' @param centroids A character vector indicating the centroids
#' @param r A numeric correlation threshold
#' @param selected_edges A character vector indicating which edges will be considered
#' @param targets  A character vector indicating the targets
#' @return A character vector with edges selected in this bootstrap
funboot<-function(seed=0,inmat,centroids,r,selected_edges,targets){
  bootsigedges<-bootmat(inmat,centroids,r,seed=seed)

  # Get surviving edges
  filtered<-bootsigedges[bootsigedges[,1]%in%centroids,]
  rm(bootsigedges)
  filtered[,1]<-as.character(filtered[,1])
  filtered[,2]<-as.character(filtered[,2])
  filtered[,3]<-as.numeric(as.character(filtered[,3]))
  rownames(filtered)<-paste0(filtered[,1],"_",filtered[,2])
  filtered<-filtered[intersect(rownames(filtered),selected_edges),]

  # Test all edges triplets for winners
  # winners<-matrix(nrow=0,ncol=3)
  # for(tg in targets){
  #   tf_candidates<-filtered[filtered[,2]==tg,]
  #   tf_candidate<-tf_candidates[which.max(abs(tf_candidates[,3])),]
  #   winners<-rbind(winners,tf_candidate)
  # }
  colnames(filtered)<-c("centroid","tg","cor")
  winners<-as.matrix(filtered %>% group_by(tg) %>% filter(abs(cor)==maxabs(cor)))
  rownames(winners)<-paste0(winners[,1],"_",winners[,2])

  # Return surviving edges in bootstrap
  return(rownames(winners))
}

#' r2p Convert Correlation Coefficient to P-value
#' @param r the correlation coefficient
#' @param n the number of samples
#' @return a numeric p-value
#' @examples
#' r2p(r=0.4,n=20) # 0.08
#' @export
r2p<-function(r,n){
  t<-(r*sqrt(n-2)) / sqrt(1-(r^2))
  p<-2*pt(t,df=n-2,lower.tail=FALSE)
  return(p)
}

#' p2r Convert a P-value to the corresponding Correlation Coefficient
#' @param p the p-value
#' @param n the number of samples
#' @return a correlation coefficient
#' @examples
#' p2r(p=0.08,n=20)
#' @export
p2r<-function(p,n){
  t<-qt(p/2,df=n-2,lower.tail=FALSE)
  r<-sqrt((t^2)/(n-2+t^2))
  return(r)
}


