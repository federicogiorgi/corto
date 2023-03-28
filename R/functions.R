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
#' @param ci logical. If TRUE, confidence intervals of linear regression are
#' shown at 95 percent confidence.
#' @param ... Arguments to be passed to the core _plot_ function (if a new plot is created)
#' @return A plot
#' @examples
#' x<-setNames(rnorm(200),paste0("var",1:200))
#' y<-setNames(rnorm(210),paste0("var",11:220))
#' scatter(x,y,xlab="Variable x",ylab="Variable y",main="Scatter plot by corto package",ci=TRUE)
#' @export
scatter<-function(x,y,method="pearson",threshold=0.01,showLine=TRUE,grid=TRUE,bgcol=FALSE,pch=20,
                  subtitle=NULL,extendXlim=FALSE,ci=FALSE,...){
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
        if(ccp<0.01){
            vv<-format(ccp,scientific=TRUE)
            v1<-gsub("e.+","",vv)
            v2<-gsub(".+e","",vv)
            v2<-gsub("-0+","-",v2)
            v2<-gsub("\\+0","+",v2)
            v2<-gsub("\\++","",v2)
            bq<-as.expression(bquote("CC="~.(cccor)~" (p="~.(v1)~x~10^.(v2)~")"))
        } else{
            bq<-mtext(paste0("CC=",cccor," (p=",ccp,")"),cex=0.7)
        }
        mtext(bq,cex=0.7)
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
    if(ci){
        lm2<-lm(y~x)
        newx=data.frame(x=seq(min(x),max(x),length.out=length(x)))
        confInterval=predict(lm2,newdata=data.frame(x=newx),interval='confidence',level=0.95)
        matlines(newx$x,confInterval[,2:3],col='dodgerblue1',lty=2)
    }
}

#' scinot - Convert a number to a scientific notation expression
#'
#' This function will convert any numeric vector
#'
#' @param v The input numeric object. It can be a single value or a vector
#' @param digits An integer indicating how many significant digits to show. Default is 3.
#' @return An object of class _expression_.
#' @examples
#' # Usage on single value
#' scinot(0.00000543)
#' # Demonstration on a vector
#' numbers<-c(3.456e-12,0.00901,5670000,-3.16e18,0.000004522,rnorm(5,sd=0.0000001))
#' plot(0,xlim=c(0,10),ylim=c(0,10),type="n")
#' text(c(2,6),c(10,10),labels=c("Before","After"),font=2)
#' for(i in 10:1){
#'     text(c(2,6),c(i-1,i-1),labels=c(numbers[i],scinot(numbers)[i]))
#' }
#' @export
scinot<-function(v,digits=3){
    v<-signif(v,digits)
    vv<-format(v,scientific=TRUE)
    v1<-gsub("e.+","",vv)
    v2<-gsub(".+e","",vv)
    v2<-gsub("-0+","-",v2)
    v2<-gsub("\\+0","+",v2)
    v2<-gsub("\\++","",v2)

    vexpr<-vector("expression",length(v))
    for(i in 1:length(vv)){
        bq<-as.expression(bquote(.(v1[i])~x~10^.(v2[i])))
        vexpr[i]<-bq
    }
    return(vexpr)
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


#' barplot2 - Bar plot with upper error bars
#' @param values A matrix of values
#' @param errors A matrix of values for upper error bar
#' @param ... Arguments to be passed to the core _barplot_ function
#' @return A plot
#' @examples
#' values<-matrix(rnorm(10*4,mean=10),nrow=4,ncol=10)
#' errors<-matrix(runif(10*4),nrow=4,ncol=10)
#' colnames(values)<-colnames(errors)<-LETTERS[1:10]
#' barplot2(values,errors,main="Bar plot with error bars")
#' @export
barplot2<-function(values,errors,...){
    if(!is.matrix(values)&!is.matrix(errors)){
        values<-t(as.matrix(values))
        errors<-t(as.matrix(errors))
    }
    sums<-values+errors
    bp<-barplot(values,beside=TRUE,ylim=c(0,1.1*max(sums)),...)
    for(i in 1:nrow(bp)){
        for(j in 1:ncol(bp)){
            pos<-bp[i,j]
            m<-values[i,j]
            s<-errors[i,j]
            arrows(pos,m+s,pos,m,angle=90,code=1,lwd=2,length=0.06)
        }
    }
}

# Function boxOverlap (needed for textrepel)
boxOverlap<-function(x1,y1,sw1,sh1,boxes){
    overlap<-TRUE
    i<-0
    while(i<length(boxes)){
        i<-i+1
        vertices<-boxes[[i]]
        x2<-vertices[1]
        y2<-vertices[2]
        sw2<-vertices[3]
        sh2<-vertices[4]
        if(x1<x2){
            overlap<-(x1+sw1)>x2
        }else{
            overlap<-(x2+sw2)>x1
        }
        if(y1<y2){
            overlap<-(overlap&&((y1+sh1)>y2))
        } else {
            overlap<-(overlap&&((y2+sh2)>y1))
        }
        if(overlap){
            return(TRUE)
        }
    }
    return(FALSE)
}


#' textrepel - Plot text with non-overlapping labels
#'
#' This function plots text with x and y coordinates, forcing overlapping labels to not overlap
#'
#' @param x A numeric vector of x coordinates
#' @param y A numeric vector of y coordinates (must have the same length of x)
#' @param labels A vector of labels associated with x and y (must have the same length of x)
#' @param padding A character object specifying left and right padding for words. Default is a
#' single whitespace " "
#' @param rstep Decimal numeric specifying the lateral step length for label distancing.
#' Default is 0.1
#' @param tstep Decimal numeric specifying the theta step length for label distancing.
#' Default is 0.1
#' @param vertical Boolean. If FALSE (default), the labels are plotted horizontally. If TRUE,
#' vertically
#' @param textSize Numeric. Size of text. Default is 1
#' @param showLines Boolean. Whether to show lines connecting displaced labels to their original
#' plot. Default is TRUE
#' @param lineColor String indicating the color of the connecting line
#' @param lineWidth Numeric indicating the width of the connecting line
#' @param showPoints Boolean. Whether to show points over original x-y coordinates
#' @param pointColor String indicating the color of the point
#' @param pointSize Numeric indicating the size of the point
#' @param pointPch Integer applying to shape of points. Default is 16 (filled circle)
#' @param add Boolean. If FALSE (default), a new plot is generated. If TRUE, the textrepel labels
#' are plotted over the existing plot
#' @param ... Arguments to be passed to the core _plot_ function
#' @return A plot
#' @examples
#' # Simple example, generating a new plot, taking care of some overlapping labels
#' set.seed(1)
#' x<-rnorm(100)
#' y<-abs(x)+rnorm(100)
#' names(x)<-names(y)<-paste0("OBJ",1:length(x))
#' labels<-names(x)
#' textrepel(x,y,labels)
#' # More advanced example, adding textrepel over an existing plot
#' set.seed(1)
#' x<-rnorm(1000)
#' y<-abs(x)+rnorm(1000)
#' names(x)<-names(y)<-paste0("GENE",1:length(x))
#' labels<-names(x)
#' plot(x,y,pch=16,col="#00000066",xlim=1.3*c(min(x),max(x)))
#' subset1<-which(x<(-2.2))
#' textrepel(x[subset1],y[subset1],labels[subset1],add=TRUE,pointCol="cornflowerblue")
#' subset2<-which(x>(+2.2))
#' textrepel(x[subset2],y[subset2],labels[subset2],add=TRUE,pointCol="salmon")


#' @export
textrepel<-function(x,y,
                    labels=NULL,
                    padding=" ",
                    rstep=0.1,
                    tstep=0.1,
                    vertical=FALSE,
                    textSize=1,
                    showLines=TRUE,
                    lineColor="#00000066",
                    lineWidth=2,
                    showPoints=TRUE,
                    pointColor="#00000033",
                    pointSize=2,
                    pointPch=16,
                    add=FALSE,
                    ...
){
    # Sanity check
    if((length(x)!=length(y))|(length(x)!=length(labels))){
        stop("Length of x, y and labels differ!")
    }

    # Limit parameters
    xLimits<-c(-Inf,Inf)
    yLimits<-c(-Inf,Inf)

    # Create new plot, if necessary
    if(!add) {
        plot(x,y,type="n",...)
    }

    # Add padding
    labels<-paste0(padding,labels,padding)

    # Define which letters are especially tall
    peculiarLetters<-"g|j|p|q|y"
    n<-length(labels)
    sdx<-sd(x,na.rm=TRUE)
    sdy<-sd(y,na.rm=TRUE)
    if(length(textSize)==1){textSize<-rep(textSize,n)}
    if(length(vertical)==1){vertical<-rep(vertical, n)}

    # The boxes loop part
    boxes<-list()
    for (i in 1:length(labels)) {
        rotLabel<-vertical[i]
        r<-0
        theta<-runif(1,0,2*pi)
        x1<-x0<-x[i]
        y1<-y0<-y[i]
        wid<-strwidth(labels[i],cex=textSize[i])
        ht<-1.2*strheight(labels[i],cex=textSize[i])
        # Increase height when peculiar letters are present
        if (grepl(peculiarLetters,labels[i])){
            ht<-ht+ht*0.2
        }
        if (rotLabel){
            tmp<-ht
            ht<-wid
            wid<-tmp
        }
        isBoxOverlapped<-TRUE
        while(isBoxOverlapped){
            if (!boxOverlap(x1-0.5*wid,y1-0.5*ht,wid,ht,boxes)&&
                (x1-0.5*wid>xLimits[1]&&y1-0.5*ht>yLimits[1])&&
                (x1+0.5*wid<xLimits[2]&&y1+0.5*ht<yLimits[2])) {
                boxes[[length(boxes)+1]]<-c(x1-0.5*wid,y1-0.5*ht,wid,ht)
                isBoxOverlapped<-FALSE
            } else{
                theta<-theta+tstep
                r<-r+rstep*tstep/(2*pi)
                x1<-x0+sdx*r*cos(theta)
                y1<-y0+sdy*r*sin(theta)
            }
        }
    }
    output<-do.call(rbind,boxes)
    colnames(output)<-c("x","y","width","ht")
    rownames(output)<-labels
    coords<-output


    # Back to plotting
    if(showLines) {
        for (i in seq_len(length(x))){
            xl<-coords[i,1]
            yl<-coords[i,2]
            w<-coords[i,3]
            h<-coords[i,4]
            if(showPoints){
                points(x[i],y[i],pch=pointPch,col=pointColor,cex=pointSize)
            }
            if(x[i]<xl||x[i]>xl+w||y[i]<yl||y[i]>yl+h){
                nx<-xl+0.5*w
                ny<-yl+0.5*h
                lines(c(x[i],nx),c(y[i],ny),col=lineColor,lwd=lineWidth)
            }
        }
    }
    for(i in 1:nrow(coords)){
        text(coords[i,1]+0.5*coords[i,3],coords[i,2]+0.5*coords[i,4],labels[i],cex=textSize[i],
             srt=vertical[i]*90,font=2)
    }
}

