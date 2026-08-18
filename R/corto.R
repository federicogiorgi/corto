#' Calculate a regulon from a data matrix
#'
#' This function applies Correlation and DPI to generate a robust
#' regulon object based on the input data matrix and the selected centroids.
#'
#' @param inmat Input matrix, with features (e.g. genes) as rows and samples
#' as columns
#' @param centroids A character vector indicating which features (e.g. genes)
#' to consider as centroids (a.k.a. Master Regulators) for DPI
#' @param nbootstraps Number of bootstraps to be performed. Default is 100
#' @param p The p-value threshold for correlation significance (by default 1E-30)
#' @param nthreads The number of threads to use for bootstrapping. Default is 1
#' @param boot_threshold The fraction of bootstraps in which the edge should appear
#'  to be included in the final network. It can be any number between 0.0 and 1.0.
#'  Default is 0.0.
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
corto<-function(inmat,centroids,nbootstraps=100,p=1E-30,nthreads=1,verbose=FALSE,cnvmat=NULL,
                boot_threshold=0.0){
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
      stop("No columns in common between cnvmat and inmat")
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
  inmat<-inmat[keep,,drop=FALSE]
  if(verbose){
    message("Removed ",nrow-length(keep)," features with zero variance")
  }
  nrow<-nrow(inmat)
  centroids<-intersect(rownames(inmat),centroids)
  targets<-setdiff(rownames(inmat),centroids)

  ## Functions required by the bootstrap jobs ----
  # Significant centroid-target correlations, on a samples x features matrix
  corhits<-function(tmat,centroids,r){
    targets<-setdiff(colnames(tmat),centroids)
    cormat<-stats::cor(tmat[,centroids,drop=FALSE],tmat[,targets,drop=FALSE])
    hits<-which(abs(cormat)>=r,arr.ind=TRUE)
    list(
      centroid=rownames(cormat)[hits[,1]],
      tg=colnames(cormat)[hits[,2]],
      cor=as.numeric(as.character(cormat[hits[,1]+nrow(cormat)*(hits[,2]-1)]))
    )
  }
  # DPI: which edges hold the strongest correlation for their target
  dpiwin<-function(tg,cor){
    a<-abs(cor)
    o<-order(tg,-a,method="radix")
    g<-tg[o]
    first<-!duplicated(g)
    mx<-rep(a[o][first],diff(c(which(first),length(g)+1L)))
    win<-logical(length(a))
    win[o]<-a[o]==mx
    win
  }
  # Function to calculate DPI on a bootstrapped matrix
  funboot<-function(seed=0,tmat,centroids,r,selected_edges){
    set.seed(seed)
    boot<-corhits(tmat[sample(nrow(tmat),replace=TRUE),,drop=FALSE],centroids,r)
    edges<-paste0(boot$centroid,"_",boot$tg)
    inset<-edges%in%selected_edges
    edges[inset][dpiwin(boot$tg[inset],boot$cor[inset])]
  }
  # Bootstraps run in their own environment, so the input matrix is shipped once
  bootenv<-new.env(parent=baseenv())
  environment(corhits)<-bootenv
  environment(dpiwin)<-bootenv
  environment(funboot)<-bootenv
  bootenv$corhits<-corhits
  bootenv$dpiwin<-dpiwin

  # Calculating pairwise correlations
  if(verbose){
    message("Calculating pairwise correlations")
  }
  tmat<-t(inmat)
  rm(inmat)
  edges<-corhits(tmat,centroids,r)
  if(length(edges$centroid)==0){
    stop("No edges passed the initial threshold, try a less stringent p")
  }

  # Extract all triplets TF-TF-TG
  if(verbose){
    message("Initial testing of triplets for DPI")
  }
  selected_edges<-paste0(edges$centroid,"_",edges$tg)
  if(verbose){
    message(length(selected_edges)," edges passed the initial threshold")
  }
  selected_nodes<-unique(c(edges$centroid,edges$tg))
  centroids<-intersect(centroids,selected_nodes)
  targets<-setdiff(selected_nodes,centroids)
  if(verbose){
    message("Building DPI network from ",length(centroids)," centroids and ",
            length(targets)," targets")
  }

  # DPI: Test all edges triplets for winners
  # Appearing in the original matrix counts as one evidence
  occurrences<-as.numeric(dpiwin(edges$tg,edges$cor))

  # Now run bootstraps to check the number of wins
  if(verbose){
    message("Running ",nbootstraps," bootstraps with ",nthreads," thread(s)")
  }

  # The pbapply snippet ----
  cl<-parallel::makeCluster(nthreads)
  winnerlist<-unlist(pblapply(cl=cl,
                              X=1:nbootstraps,
                              FUN=funboot,
                              tmat=tmat,centroids=centroids,r=r,selected_edges=selected_edges
  ))
  parallel::stopCluster(cl)

  # Add occurrences
  add<-table(winnerlist)
  hitrows<-match(names(add),selected_edges)
  occurrences[hitrows]<-occurrences[hitrows]+as.numeric(add)

  # Likelihood based on bootstrap occurrence
  if(verbose){
    message("Calculating edge likelihood")
  }
  likelihood<-occurrences/(nbootstraps+1)
  keep<-likelihood>boot_threshold

  # Generate regulon object
  if(verbose){
    message("Generating regulon object")
  }
  cens<-edges$centroid[keep]
  tgs<-edges$tg[keep]
  cors<-edges$cor[keep]
  liks<-likelihood[keep]
  rows<-split(seq_along(cens),cens)
  regulon<-list()
  for(tf in unique(cens)){
    i<-rows[[tf]]
    regulon[[tf]]<-list(
      tfmode=setNames(cors[i],tgs[i]),
      likelihood=liks[i]
    )
  }

  class(regulon)<-"regulon"
  return(regulon)
}



