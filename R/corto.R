#' Calculate a regulon from a data matrix
#'
#' This function applies Spearman Correlation and DPI to generate a robust
#' regulon object based on the input data matrix and the selected centroids.
#'
#' @param inmat Input matrix, with features (e.g. genes) as rows and samples
#' as columns
#' @param centroids A character vector indicating which features (e.g. genes)
#' to consider as centroids (a.k.a. Master Regulators) for DPI
#' @param nbootstraps Number of bootstraps to be performed. Default is 100
#' @param p The p-value threshold for correlation significance (by default 1E-30)
#' @return A list (object of class regulon), where each element is a centroid
#' \itemize{
#'   \item tfmode: a named vector containing correlation coefficients between
#'   features and the centroid
#'   \item likelihood: a numeric vector indicating the likelihood of interaction
#' }
#' @examples
#' # Load data matrix inmat (from TCGA mesothelioma project)
#' load(system.file("extdata","inmat.rda",package="corto",mustWork=TRUE))
#' # Load centroids
#' load(system.file("extdata","centroids.rda",package="corto",mustWork=TRUE))
#' # Run corto
#' regulon <- corto(inmat,centroids=centroids)
#' @export
corto<-function(inmat,centroids,nbootstraps=100,p=1E-30){
  # Analytical inference of threshold
  ncol<-ncol(inmat)
  nrow<-nrow(inmat)
  message("Input Matrix has ",ncol," samples and ",nrow," features")
  r<-p2r(p=p,n=ncol)
  message("Correlation Coefficient Threshold is: ",r)

  # Filtering zero variance
  allvars<-apply(inmat,1,var)
  keep<-names(allvars)[allvars>0]
  inmat<-inmat[keep,]
  message("Removed ",nrow-length(keep)," features with zero variance")
  nrow<-nrow(inmat)
  centroids<-intersect(rownames(inmat),centroids)
  targets<-setdiff(rownames(inmat),centroids)

  # Calculating pairwise correlations
  message("Calculating pairwise correlations")
  sigedges<-fcor(inmat,centroids,r)

  # Extract all triplets TF-TF-TG
  message("Initial testing of triplets for DPI")
  # Remove edges which have no TF
  filtered<-sigedges[sigedges[,1]%in%centroids,]
  rm(sigedges)
  filtered[,1]<-as.character(filtered[,1])
  filtered[,2]<-as.character(filtered[,2])
  filtered[,3]<-as.numeric(as.character(filtered[,3]))
  rownames(filtered)<-paste0(filtered[,1],"_",filtered[,2])
  selected_edges<-rownames(filtered)
  message(length(selected_edges)," edges passed the initial threshold")
  selected_nodes<-unique(c(filtered[,1],filtered[,2]))
  centroids<-intersect(centroids,selected_nodes)
  targets<-setdiff(selected_nodes,centroids)
  message("Building network with ",length(centroids)," centroids and ",
          length(targets)," targets")

  # DPI: Test all edges triplets for winners
  winners<-matrix(nrow=0,ncol=3)
  for(tg in targets){
    tf_candidates<-filtered[filtered[,2]==tg,]
    tf_candidate<-tf_candidates[which.max(abs(tf_candidates[,3])),]
    winners<-rbind(winners,tf_candidate)
  }
  occ<-cbind(filtered,rep(0,nrow(filtered)))
  #occ[rownames(winners),4]<-100
  rm(winners,filtered)
  colnames(occ)[4]<-"occurrences"

  # Now run bootstraps to check the number of wins
  message("Running ",nbootstraps," bootstraps")
  pb<-txtProgressBar(0,nbootstraps,style=3)
  for(i in 1:nbootstraps){
    bootsigedges<-bootfun(inmat,centroids,r,seed=i)

    # Get surviving edges
    filtered<-bootsigedges[bootsigedges[,1]%in%centroids,]
    rm(bootsigedges)
    filtered[,1]<-as.character(filtered[,1])
    filtered[,2]<-as.character(filtered[,2])
    filtered[,3]<-as.numeric(as.character(filtered[,3]))
    rownames(filtered)<-paste0(filtered[,1],"_",filtered[,2])
    filtered<-filtered[intersect(rownames(filtered),selected_edges),]

    # Test all edges triplets for winners
    winners<-matrix(nrow=0,ncol=3)
    for(tg in targets){
      tf_candidates<-filtered[filtered[,2]==tg,]
      tf_candidate<-tf_candidates[which.max(abs(tf_candidates[,3])),]
      winners<-rbind(winners,tf_candidate)
    }
    # Update original count for occurrences
    occ[rownames(winners),4]<-occ[rownames(winners),4]+1

    # Update progress bar
    setTxtProgressBar(pb,i)
  }

  # Likelihood based on bootstrap occurrence
  message("Calculating edge likelihood")
  occ$likelihood<-occ$occurrences/nbootstraps

  # Generate regulon object
  message("Generating regulon object")
  regulon<-list()
  for(tf in unique(occ[,1])){
    targets<-occ[occ[,1]==tf,2]
    corr<-occ[occ[,1]==tf,3]
    likelihood<-occ[occ[,1]==tf,5]
    tfmode<-setNames(corr,targets)
    regulon[[tf]]<-list(
      tfmode=tfmode,
      likelihood=likelihood
    )
  }

  class(regulon)<-"regulon"
  return(regulon)
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
  p<-2*pt(t,df=n-2,lower=FALSE)
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
  t<-qt(p/2,df=n-2,lower=FALSE)
  r<-sqrt((t^2)/(n-2+t^2))
  return(r)
}

