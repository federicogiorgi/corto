#' Apply the Data Processing Inequality to a set of centroid-target edges
#'
#' Decides which centroid-target edges survive the Data Processing Inequality
#' (DPI). Used by \code{corto()} on the original data and in every bootstrap, and
#' usable on its own with any list of correlation edges.
#'
#' @param centroid Character vector, the centroid (regulator) of each edge
#' @param target Character vector, the target of each edge
#' @param cor Numeric vector, the correlation of each edge
#' @param rule \code{"target"} (default, the only rule of \code{corto()} up to
#'  version 1.3.1): for each target, only the centroid(s) with the highest absolute
#'  correlation survive. \code{"triangle"}: the ARACNe rule; the edge between
#'  centroid \emph{a} and target \emph{c} is removed only if another centroid
#'  \emph{b} closes a triangle in which \emph{a}-\emph{c} is the weakest edge,
#'  i.e. if both |cor(b,c)| and |cor(a,b)| are larger than |cor(a,c)|. Centroids
#'  that are not correlated with each other keep their edges to a shared target.
#' @param centroid_cor Correlation matrix between centroids, with row and column
#'  names; required for \code{rule="triangle"}. Only absolute values are used.
#' @return A logical vector, \code{TRUE} for the edges that survive
#' @examples
#' # Two centroids correlated with the same target. Under "target" only the
#' # stronger survives; under "triangle" both survive unless they are more
#' # correlated with each other than the weaker one is with the target.
#' cc <- matrix(c(1, 0.1, 0.1, 1), 2, dimnames = list(c("A", "B"), c("A", "B")))
#' dpi_edges(c("A", "B"), c("g", "g"), c(0.6, 0.4), rule = "target")
#' dpi_edges(c("A", "B"), c("g", "g"), c(0.6, 0.4), rule = "triangle", centroid_cor = cc)
#' @export
dpi_edges<-function(centroid,target,cor,rule=c("target","triangle"),centroid_cor=NULL){
  rule<-match.arg(rule)
  a<-abs(cor)
  if(rule=="target"){
    # Winner per target: radix sort by target and decreasing |cor|
    o<-order(target,-a,method="radix")
    g<-target[o]
    first<-!duplicated(g)
    mx<-rep(a[o][first],diff(c(which(first),length(g)+1L)))
    win<-logical(length(a))
    win[o]<-a[o]==mx
    return(win)
  }
  if(is.null(centroid_cor)){
    stop("rule=\"triangle\" needs centroid_cor, the correlation matrix between centroids")
  }
  centroid_cor<-abs(centroid_cor)
  keep<-rep(TRUE,length(a))
  for(i in split(seq_along(a),target)){
    if(length(i)<2) next
    ai<-a[i]
    # [x,y] is TRUE when centroid y beats x on this target and is more
    # correlated with x than x is with the target
    beaten<-(centroid_cor[centroid[i],centroid[i],drop=FALSE]>ai)&outer(ai,ai,"<")
    keep[i]<-rowSums(beaten)==0
  }
  keep
}

#' Calculate a regulon from a data matrix
#'
#' This function applies Correlation and DPI to generate a robust
#' regulon object based on the input data matrix and the selected centroids.
#'
#' @details
#' All centroid-target Pearson correlations are calculated on \code{inmat} as
#' provided (no rank transformation is applied: rank the matrix beforehand to
#' work on Spearman correlations), and the edges whose absolute correlation
#' passes the threshold derived from \code{p} are kept. The Data Processing
#' Inequality (DPI) is then applied to these edges on the original matrix and on
#' every bootstrapped matrix, and the likelihood of an edge is the number of
#' times it survives DPI divided by \code{nbootstraps+1}.
#'
#' Two DPI rules are available (see \code{\link{dpi_edges}}):
#' \itemize{
#'   \item \code{"target"} (default): for each target, only the edge from the
#'   centroid with the highest absolute correlation survives. This assumes one
#'   regulator per target: co-regulators of the same target are only retained
#'   through bootstrap variability.
#'   \item \code{"triangle"}: the ARACNe rule. An edge is removed only if it is
#'   the weakest edge of a centroid-centroid-target triangle, so co-regulators
#'   that are not correlated with each other keep their edges. The centroid-centroid
#'   correlations are calculated on the same (bootstrapped) matrix. Networks are
#'   larger and closer to those of ARACNe, at a considerably higher computing time.
#' }
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
#' @param dpi The DPI rule, passed to \code{\link{dpi_edges}}: \code{"target"}
#'  (default, one winning centroid per target, as in previous versions) or
#'  \code{"triangle"} (the ARACNe rule over centroid-centroid-target triangles).
#'  See Details.
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
                boot_threshold=0.0,dpi=c("target","triangle")){
  dpi<-match.arg(dpi)
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
    cnvmat<-cnvmat[commonrows,commoncols,drop=FALSE]
    inmat<-inmat[commonrows,commoncols,drop=FALSE]

    # Correct inmat based on cnvmat (only for targets)
    targetmat<-inmat[setdiff(rownames(cnvmat),centroids),,drop=FALSE]
    cnvtargetmat<-cnvmat[setdiff(rownames(cnvmat),centroids),,drop=FALSE]
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
  # DPI on a set of edges, the centroid-centroid correlations on the same
  # samples x features matrix are only needed by the "triangle" rule
  dpiedges<-dpi_edges
  dpimat<-function(edges,tmat,dpi){
    if(length(edges$cor)==0){
      return(logical(0))
    }
    cc<-NULL
    if(dpi=="triangle"){
      cc<-stats::cor(tmat[,unique(edges$centroid),drop=FALSE])
    }
    dpiedges(edges$centroid,edges$tg,edges$cor,rule=dpi,centroid_cor=cc)
  }
  # Function to calculate DPI on a bootstrapped matrix
  funboot<-function(seed=0,tmat,centroids,r,selected_edges,dpi){
    set.seed(seed)
    bmat<-tmat[sample(nrow(tmat),replace=TRUE),,drop=FALSE]
    boot<-corhits(bmat,centroids,r)
    edges<-paste0(boot$centroid,"_",boot$tg)
    inset<-edges%in%selected_edges
    sub<-lapply(boot,function(x) x[inset])
    edges[inset][dpimat(sub,bmat,dpi)]
  }
  # Bootstraps run in their own environment, so the input matrix is shipped once
  bootenv<-new.env(parent=baseenv())
  environment(corhits)<-bootenv
  environment(dpiedges)<-bootenv
  environment(dpimat)<-bootenv
  environment(funboot)<-bootenv
  bootenv$corhits<-corhits
  bootenv$dpiedges<-dpiedges
  bootenv$dpimat<-dpimat

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

  # Edges entering the DPI step
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
    if(dpi=="triangle"){
      message("Using the triangle DPI rule")
    }
  }

  # DPI on the original matrix
  # Appearing in the original matrix counts as one evidence
  occurrences<-as.numeric(dpimat(edges,tmat,dpi))

  # Now run bootstraps to check the number of wins
  if(verbose){
    message("Running ",nbootstraps," bootstraps with ",nthreads," thread(s)")
  }

  # The pbapply snippet ----
  cl<-parallel::makeCluster(nthreads)
  winnerlist<-unlist(pblapply(cl=cl,
                              X=1:nbootstraps,
                              FUN=funboot,
                              tmat=tmat,centroids=centroids,r=r,selected_edges=selected_edges,
                              dpi=dpi
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



