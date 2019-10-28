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
#' @param nthreads The number of threads to use for bootstrapping. Default is 1
#' @param verbose Logical. Whether to print progress messages. Default is FALSE
#' @param cnvmat An optional matrix with copy-number variation data. If specified, the program
#' will calculate linear regression between the gene expression data in the input matrix (exp)
#' and the cnv data, and target profiles will be transformed to the residuals
#' of each linear model exp~cnv. Default is NULL
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
#' regulon <- corto(inmat,centroids=centroids,nthreads=2,nbootstraps=10,verbose=TRUE)
#'
#' # In a second example, a CNV matrix is provided. The analysis will be run only
#' # for the features (rows) and samples (columns) present in both matrices
#' load(system.file("extdata","cnvmat.rda",package="corto",mustWork=TRUE))
#' regulon <- corto(inmat,centroids=centroids,nthreads=2,nbootstraps=6,verbose=TRUE,cnvmat=cnvmat,
#' p=1e-8)
#' @export
corto<-function(inmat,centroids,nbootstraps=100,p=1E-30,nthreads=1,verbose=FALSE,cnvmat=NULL){
  if(sum(is.na(inmat))>0){
    stop("Input matrix contains NA fields")
  }
  if(sum(is.na(rownames(inmat)))>0){
    stop("Row names of inmat contain NA values")
  }
  if(any(abs(inmat)==Inf)){
    stop("Input matrix contains Inf fields")
  }
  # Initial check for the presence of CNV data
  if(!is.null(cnvmat)){
    if(verbose){
      message("Correcting input data with user-provided CNV data...")
    }
    if(sum(is.na(cnvmat))>0){
      stop("Input cnvmat contains NA fields")
    }
    # Intersect the two matrices
    commonrows<-intersect(rownames(cnvmat),rownames(inmat))
    if(length(commonrows)<=1){
      stop("One or less rows in common between cnvmat and inmat")
    }
    commoncols<-intersect(colnames(cnvmat),colnames(inmat))
    if(length(commoncols)==0){
      stop("One or less columns in common between cnvmat and inmat")
    }
    cnvmat<-cnvmat[commonrows,commoncols]
    inmat<-inmat[commonrows,commoncols]

    # Correct inmat based on cnvmat (only for targets)
    targetmat<-inmat[setdiff(rownames(cnvmat),centroids),]
    cnvtargetmat<-cnvmat[setdiff(rownames(cnvmat),centroids),]
    if(verbose){
      message("Applying residual calculation for ",ncol(targetmat)," samples and ",nrow(targetmat)," target features")
    }
    corrected<-t(apply(cbind(targetmat,cnvtargetmat),1,function(inline){
      y<-inline[1:ncol(targetmat)]
      x<-inline[(ncol(targetmat)+1):(length(inline))]
      lm1<-lm(y~x)
      return(lm1$residuals)
    }))
    inmat[rownames(targetmat),]<-corrected[rownames(targetmat),]
  }

  # Analytical inference of threshold
  ncol<-ncol(inmat)
  nrow<-nrow(inmat)
  if(verbose){
    message("Input Matrix has ",ncol," samples and ",nrow," features")
  }
  r<-p2r(p=p,n=ncol)
  if(verbose){
    message("Correlation Coefficient Threshold is: ",r)
  }
  # Filtering zero variance
  allvars<-apply(inmat,1,var)
  keep<-names(allvars)[allvars>0]
  inmat<-inmat[keep,]
  if(verbose){
    message("Removed ",nrow-length(keep)," features with zero variance")
  }
  nrow<-nrow(inmat)
  centroids<-intersect(rownames(inmat),centroids)
  targets<-setdiff(rownames(inmat),centroids)

  # Calculating pairwise correlations
  if(verbose){
    message("Calculating pairwise correlations")
  }
  sigedges<-fcor(inmat,centroids,r)

  # Extract all triplets TF-TF-TG
  if(verbose){
    message("Initial testing of triplets for DPI")
  }
  # Remove edges which have no TF
  filtered<-sigedges[sigedges[,1]%in%centroids,]
  rm(sigedges)
  filtered[,1]<-as.character(filtered[,1])
  filtered[,2]<-as.character(filtered[,2])
  filtered[,3]<-as.numeric(as.character(filtered[,3]))
  rownames(filtered)<-paste0(filtered[,1],"_",filtered[,2])
  selected_edges<-rownames(filtered)
  if(verbose){
    message(nrow(filtered)," edges passed the initial threshold")
  }
  selected_nodes<-unique(c(filtered[,1],filtered[,2]))
  centroids<-intersect(centroids,selected_nodes)
  targets<-setdiff(selected_nodes,centroids)
  if(verbose){
    message("Building DPI network from ",length(centroids)," centroids and ",
            length(targets)," targets")
  }

  # DPI: Test all edges triplets for winners
  # winners<-matrix(nrow=0,ncol=3)
  # for(tg in targets){ # TO DO an apply
  #   tf_candidates<-filtered[filtered[,2]==tg,]
  #   tf_candidate<-tf_candidates[which.max(abs(tf_candidates[,3])),]
  #   winners<-rbind(winners,tf_candidate)
  # }
  colnames(filtered)<-c("centroid","tg","cor")
  tg<-filtered[,"tg"]
  winners<-as.matrix(filtered %>% group_by(tg) %>% filter(abs(cor)==maxabs(cor)))
  rownames(winners)<-paste0(winners[,1],"_",winners[,2])
  occ<-cbind(filtered,rep(0,nrow(filtered)))
  # Appearing in the original matrix counts as one evidence
  occ[rownames(winners),4]<-1
  rm(winners,filtered)
  colnames(occ)[4]<-"occurrences"

  # Now run bootstraps to check the number of wins
  if(verbose){
    message("Running ",nbootstraps," bootstraps with ",nthreads," thread(s)")
  }


  ## Run the bootstraps in multithreading
  # Functions required by the bootstrap jobs ----
  # Function to calculate DPI
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
    colnames(filtered)<-c("centroid","tg","cor")
    tg<-filtered[,"tg"]
    winners<-as.matrix(filtered %>% group_by(tg) %>% filter(abs(cor)==maxabs(cor)))
    rownames(winners)<-paste0(winners[,1],"_",winners[,2])
    # Return surviving edges in bootstrap
    return(rownames(winners))
  }
  # Function to bootstrap Matrix
  bootmat<-function(inmat,centroids,r,seed=NULL){
    set.seed(seed)
    bootmat<-inmat[,sample(colnames(inmat),replace=TRUE)]
    # Calculate correlations in the bootstrapped matrix
    bootsigedges<-fcor(bootmat,centroids,r)
    return(bootsigedges)
  }
  # Function to calculate correlation in bootstraps
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
  # Max of absolute value
  maxabs<-function(x){
    max(abs(x))
  }



  # The pbapply snippet ----
  cl<-parallel::makeCluster(nthreads)
  #invisible(parallel::clusterExport(cl=cl,varlist=c("nthreads")))
  invisible(parallel::clusterEvalQ(cl= cl,library(dplyr)))
  winnerlist<-unlist(pblapply(cl=cl,
                              X=1:nbootstraps,
                              FUN=funboot,
                              inmat=inmat,centroids=centroids,r=r,selected_edges=selected_edges,targets=targets
  ))
  parallel::stopCluster(cl)

  # Add occurrences
  add<-table(winnerlist)
  occ[names(add),"occurrences"]<-occ[names(add),"occurrences"]+add

  # Likelihood based on bootstrap occurrence
  if(verbose){
    message("Calculating edge likelihood")
  }
  occ$likelihood<-occ$occurrences/(nbootstraps+1)
  occ<-occ[occ$likelihood>0,]

  # Generate regulon object
  if(verbose){
    message("Generating regulon object")
  }
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



