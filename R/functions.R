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

# #' Function to bootstrap matrix
# #' @param inmat An input matrix with features as rows and samples as columns
# #' @param centroids A character vector indicating the centroids
# #' @param r A numeric correlation threshold
# #' @param seed An integer to set the random seed
# #' @return A matrix describing which edges were significant in the bootstrapped
# #' matrix according to the r correlation threshold provided
# bootmat<-function(inmat,centroids,r,seed=NULL){
#   set.seed(seed)
#   bootmat<-inmat[,sample(colnames(inmat),replace=TRUE)]
#   # Calculate correlations in the bootstrapped matrix
#   options(warn=-1)
#   bootsigedges<-fcor(bootmat,centroids,r)
#   options(warn=0)
#   return(bootsigedges)
# }

# #' Function to process bootstraps into edges
# #' @param seed An integer to set the random seed
# #' @param inmat An input matrix with features as rows and samples as columns
# #' @param centroids A character vector indicating the centroids
# #' @param r A numeric correlation threshold
# #' @param selected_edges A character vector indicating which edges will be considered
# #' @param targets  A character vector indicating the targets
# #' @return A character vector with edges selected in this bootstrap
# funboot<-function(seed=0,inmat,centroids,r,selected_edges,targets){
#   bootsigedges<-bootmat(inmat,centroids,r,seed=seed)
#
#   # Get surviving edges
#   filtered<-bootsigedges[bootsigedges[,1]%in%centroids,]
#   rm(bootsigedges)
#   filtered[,1]<-as.character(filtered[,1])
#   filtered[,2]<-as.character(filtered[,2])
#   filtered[,3]<-as.numeric(as.character(filtered[,3]))
#   rownames(filtered)<-paste0(filtered[,1],"_",filtered[,2])
#   filtered<-filtered[intersect(rownames(filtered),selected_edges),]
#
#   # Test all edges triplets for winners
#   # winners<-matrix(nrow=0,ncol=3)
#   # for(tg in targets){
#   #   tf_candidates<-filtered[filtered[,2]==tg,]
#   #   tf_candidate<-tf_candidates[which.max(abs(tf_candidates[,3])),]
#   #   winners<-rbind(winners,tf_candidate)
#   # }
#   colnames(filtered)<-c("centroid","tg","cor")
#   tg<-filtered[,"tg"]
#   winners<-as.matrix(filtered %>% group_by(tg) %>% filter(abs(cor)==maxabs(cor)))
#   rownames(winners)<-paste0(winners[,1],"_",winners[,2])
#
#   # Return surviving edges in bootstrap
#   return(rownames(winners))
# }

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


#' kmgformat - Nice Formatting of Numbers
#'
#' This function will convert thousand numbers to K, millions to M, billions
#' to G, trillions to T, quadrillions to P
#'
#' @param input A vector of values
#' @param roundParam How many decimal digits you want
#' @return A character vector of formatted numebr names
#' @examples
#' # Thousands
#' set.seed(1)
#' a<-runif(1000,0,1e4)
#' plot(a,yaxt='n')
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#'
#' # Millions to Billions
#' set.seed(1)
#' a<-runif(1000,0,1e9)
#' plot(a,yaxt='n',pch=20,col="black")
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#' @export
kmgformat <- function(input, roundParam = 1) {
  signs <- sign(input)
  signs[signs == 1] <- ""
  signs[signs == -1] <- "-"
  signs[signs == 0] <- ""
  absinput <- abs(input)
  output <- c()
  for (i in absinput) {
    if (i < 1000) {
      output <- c(output, i)
    } else if (i < 1e+06) {
      i <- round(i/1000, roundParam)
      i <- paste0(i, "K")
      output <- c(output, i)
    } else if (i < 1e+09) {
      i <- round(i/1e+06, roundParam)
      i <- paste0(i, "M")
      output <- c(output, i)
    } else if (i < 1e+12) {
      i <- round(i/1e+09, roundParam)
      i <- paste0(i, "G")
      output <- c(output, i)
    } else if (i < 1e+15) {
      i <- round(i/1e+12, roundParam)
      i <- paste0(i, "T")
      output <- c(output, i)
    } else if (i < 1e+18) {
      i <- round(i/1e+15, roundParam)
      i <- paste0(i, "P")
      output <- c(output, i)
    } else {
      output <- c(output, i)
    }
  }
  output <- paste0(signs, output)
  return(output)
}

#' scatter - XY scatter plot with extra information
#'
#' This function will plot two variables (based on their common names), calculate their
#' Coefficient of Correlation (CC), plot a linear regression line and color the background
#' if the correlation is positive (red), negative
#' (blue) or non-significant (white)
#'
#' @param x The first named vector
#' @param y The second named vector
#' @param method a character string indicating which correlation coefficient is to
#' be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param threshold a numeric value indicating the significance threshold (p-value) of the correlation,
#' in order to show a colored background. Default is 0.01.
#' @param showLine a boolean indicating if a linear regression line should be plotted. Default is
#' TRUE
#' @param grid a boolean indicating whether to show a plot grid. Default is TRUE
#' @param pch the _pch_ parameter indicating the points shape. Default is 20
#' @param subtitle NULL by default, in which case the function will print as a subtitle the correlation
#' coefficient (CC) and its pvalue. Otherwise, a user-provided string, bypassing the predefined subtitle
#' @param extendXlim logical. If TRUE, the x-axis limits are extended by a fraction (useful for
#' labeling points on the margins of the plot area). Default is FALSE
#' @param bgcol Boolean. Should a background coloring associated to significance and sign of
#' correlation be used? Default is TRUE, and it will color the background in red if the correlation
#' coefficient is positive, in blue if negative, in white if not significant (accordin to the
#' _threshold_ parameter)
#' @param ... Arguments to be passed to the core _plot_ function
#' @return A plot
#' @examples
#' x<-setNames(rnorm(200),paste0("var",1:200))
#' y<-setNames(rnorm(210),paste0("var",11:220))
#' scatter(x,y,xlab="Variable x",ylab="Variable y",main="Scatter plot by corto package")
#' @export
scatter<-function(x,y,method="pearson",threshold=0.01,showLine=TRUE,grid=TRUE,bgcol=FALSE,pch=20,
                  subtitle=NULL,extendXlim=FALSE,...){
  common<-intersect(names(x),names(y))
  x<-x[common]
  y<-y[common]
  if(!extendXlim){
    plot(x,y,pch=pch,...)
  }else{
    plot(x,y,pch=pch,xlim=1.2*c(min(x),max(x)),...)
  }
  cc<-cor.test(x,y,method=method)
  ccp<-signif(cc$p.value,3)
  cccor<-signif(cc$estimate,3)
  if(is.null(subtitle)){
    mtext(paste0("CC=",cccor," (p=",ccp,")"),cex=0.7)
  } else {
    mtext(subtitle,cex=0.7)
  }
  if(bgcol){
    if(cccor>=0){bgcol<-"#FF000033"}else{bgcol<-"#0000FF33"}
    if(ccp>threshold){bgcol<-"#FFFFFF00"}
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],col=bgcol)
  }
  if(grid){
      grid(col="gray10")
  }
  if(showLine){
    lm1<-lm(y~x)
    abline(lm1$coef)
  }
}


### Arena: Advanced Rank ENrichment Analysis
arena<-function(
  signatures,
  groups,
  sweights=NULL,
  gweights=NULL,
  minsize=1
){
  ### Convert to list if groups are not a list
  if(!is.list(groups)){
    groups<-apply(groups,1,function(x){
      hits<-colnames(groups)[x==1]
      return(hits)
    })
  }

  ### Remove from groups what is not in signature matrix
  features<-rownames(signatures)
  groups<-lapply(groups,function(x){
    y<-intersect(x,features)
    return(y)
  })

  ### Remove small groups
  groups<-groups[sapply(groups,length)>=minsize]

  ### Treat single "signature"
  if (is.null(nrow(signatures))){
    signatures <- matrix(signatures, length(signatures), 1, dimnames=list(names(signatures), "sample1"))
  }


  ### Generate dummy signature weights
  if(is.null(sweights)){
    sweights<-matrix(1,nrow=nrow(signatures),ncol=ncol(signatures))
    dimnames(sweights)<-dimnames(signatures)
  }
  if(!identical(dim(signatures), dim(sweights))){
    stop("Signatures and Signature weights must be matrices of identical size")
  }

  ### Generate dummy group weights
  if(is.null(gweights)){
    gweights<-relist(rep(1,sum(sapply(groups,length))),skeleton=groups)
  }

  ### Apply weights to group belonging
  wgroups<-gweights
  for(i in 1:length(wgroups)){
    names(wgroups[[i]])<-groups[[i]]
  }
  rm(gweights)

  ### Rank-transform columns
  ranks<-apply(signatures,2,rank,na.last="keep")
  ### Assign a 0 to signature weights where the signature was NA
  sweights[is.na(ranks)]<-0
  ### 0-1 bound ranks
  boundranks<-t(t(ranks)/(colSums(!is.na(signatures))+1))
  ### Treat bound ranks as quantiles in a gaussian distribution (0=-Inf, 1=+Inf)
  gaussian <- qnorm(boundranks)
  ### Deal with NAs
  gaussian[is.na(gaussian)]<-0

  ### Apply signature weights to the normalized distribution
  gaussian<-gaussian*sweights


  ### Next, we see how each of the groups are behaving in these normalized signatures
  ### Create a boolean matrix with ngroup columns and signaturelength rows, indicating the matches
  matches <- sapply(wgroups, function(group, allElements) {
    hereMatches<-as.integer(allElements%in%names(group))
    names(hereMatches)<-allElements
    # Weigth by group belonging
    weightedMatches<-hereMatches
    weightedMatches[names(group)]<-weightedMatches[names(group)]*group
    return(weightedMatches)
  }, allElements=rownames(gaussian))
  # And then transpose it
  matches<-t(matches)
  colnames(matches)<-rownames(signatures)

  # Number of matches per group
  groupmatches <- rowSums(matches)

  # Relative part of the signature that matches
  relativematches<-matches/groupmatches

  # This trick will overweight massively small groups with all their components highly-ranked.
  # Extreme case is with a group with one gene at the top

  # The core linear algebra operation. The true magic of rea
  enrichmentScore <- relativematches %*% gaussian

  # Finally, every enrichment is square-rooted to respect the criterion of normality
  normalizedEnrichmentScore<-enrichmentScore*sqrt(groupmatches)

  # Return output
  return(normalizedEnrichmentScore)
}

#' val2col - Convert a numeric vector into colors
#' @param z a vector of numbers
#' @param col1 a color name for the min value, default 'navy'
#' @param col2 a color name for the middle value, default 'white'
#' @param col3 a color name for the max value, default 'red3'
#' @param nbreaks Number of colors to be generated. Default is 30.
#' @param center boolean, should the data be centered? Default is TRUE
#' @param rank boolean, should the data be ranked? Default is FALSE
#' @return a vector of colors
#' @examples
#' a<-rnorm(1000)
#' cols<-val2col(a)
#' plot(a,col=cols,pch=16)
#' @export
val2col <- function(z, col1 = "navy", col2 = "white",
                    col3 = "red3", nbreaks = 1000, center = TRUE,
                    rank = FALSE) {
  isMatrix <- FALSE
  if (is.matrix(z)) {
    isMatrix <- TRUE
    oriz <- z
  }
  if (is.character(z)) {
    z <- as.numeric(as.factor(z))
  }
  if (rank) {
    z <- rank(z)
  }
  if (center) {
    extreme = round(max(abs(z)))
    breaks <- seq(-extreme, extreme, length = nbreaks)
    z <- z - mean(z)
  } else {
    breaks <- seq(min(z), max(z), length = nbreaks)
  }
  ncol <- length(breaks) - 1
  col <- gplots::colorpanel(ncol, col1, col2, col3)
  CUT <- cut(z, breaks = breaks)
  # assign colors to heights for each point
  colorlevels <- col[match(CUT, levels(CUT))]
  names(colorlevels) <- names(z)
  if (isMatrix) {
    colormatrix <- matrix(colorlevels, ncol = ncol(oriz), nrow = nrow(oriz))
    dimnames(colormatrix) <- dimnames(oriz)
    return(colormatrix)
  }
  return(colorlevels)
}
