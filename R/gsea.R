#' GSEA
#'
#' This function performs Gene Set Enrichment Analysis
#' @param reflist named vector of reference scores
#' @param set element set
#' @param method one of 'permutation' or 'pareto'
#' @param w exponent used to raise the supplied scores. Default is 1 (original
#' scores unchanged)
#' @param np Number of permutations (Default: 1000)
#' @param gsea_null a GSEA null distribution (Optional)
#' @return A GSEA object. Basically a list of s components:
#' \describe{
#' \item{ES}{The enrichment score}
#' \item{NES}{The normalized enrichment socre}
#' \item{ledge}{The items in the leading edge}
#' \item{p.value}{The permutation-based p-value}
#' }
#' @examples
#' reflist<-setNames(-sort(rnorm(1000)),paste0('gene',1:1000))
#' set<-paste0('gene',sample(1:200,50))
#' obj<-gsea(reflist,set,method='pareto',np=1000)
#' obj$p.value
#' @export
gsea <- function(reflist,
                 set,
                 method=c("permutation","pareto"),
                 np=1000,
                 w=1,
                 gsea_null=NULL) {


    # Get elements in set that are in the ref list
    set <- intersect(names(reflist), set)

    # Sort the reference list
    # Get the list order, from higher (1)to smaller (n)
    ix <- order(reflist, decreasing=TRUE)
    reflist <- reflist[ix] # Reorder the reference list

    # Initialize variables for running sum
    es <- 0
    nes <- 0
    p.value <- 1

    # Identify indexes of set within the sorted reference list
    inSet <- rep(0, length(reflist))
    inSet[which(names(reflist) %in% set)] <- 1

    ### Compute Enrichment Score
    # Compute running sum for hits
    hits<-abs(reflist*inSet) # Get the values for the elements in the set
    hits<-hits^w # Raise this score to the power of w
    score_hit <- cumsum(hits) # Cumulative sum of hits' scores
    # The cumulative sum is divided by the final  sum value
    score_hit <- score_hit / score_hit[length(score_hit)]

    # Compute running sum for non-hits
    score_miss <- cumsum(1-inSet)
    score_miss <- score_miss/score_miss[length(score_miss)]

    # The Running Score is the difference between the two scores! Hits - nonhits
    running_score <- score_hit - score_miss

    # Safety measure, in the case the random genes have all a weight of 0
    if(all(is.na(running_score))){
        running_score<-rep(0,length(running_score))
    }

    # The ES is actually the minimum or maximum Running Scores
    if(abs(max(running_score))>abs(min(running_score))){
        es<-max(running_score)
    } else {
        es<-min(running_score)
    }


    ### Identify leading edge
    # Create a vector of 0s long as the reference list
    ledge_indeces <- rep(0, length(running_score))
    # Case 1: negative ES
    if (es<0){
        peak <- which(running_score==min(running_score))[1]
        # Leading edge is stuff AFTER the peak point (ES is negative)
        ledge_indeces[peak:length(ledge_indeces)] <- 1
        ledge_indeces <- which(ledge_indeces == 1)
        ledge_names <- names(reflist[ledge_indeces])
    } else{ # Case 2: positive ES
        peak <- which(running_score==max(running_score)) # Define the peak point
        # Leading edge is stuff BEFORE the peak point (ES is positive)
        ledge_indeces[1:peak] <- 1
        ledge_indeces <- which(ledge_indeces == 1)
        ledge_names <- names(reflist[ledge_indeces])
    }





    ### Compute p-value by permutation
    if(is.null(gsea_null)){
        null_es<-null_gsea(set=set,reflist=reflist,np=np,w=w)
    } else{
        ### If a null list is provided, use it
        if(class(gsea_null)=="gsea_nullist"){
            null_es<-gsea_null[as.character(length(set))][[1]]
        }else{
            null_es<-gsea_null
        }
    }
    if (es<0){
        p.value <- sum(null_es<=es)/length(null_es)
    } else {
        p.value <- sum(null_es>=es)/length(null_es)
    }
    #    }


    # If we are in the tail, the p-value can be calculated in two ways
    if(is.na(p.value) || p.value<0.05) {
        if(p.value==0){
            p.value <- 1/np
        }
        if (method=="pareto"){
            # Extract the absolute null ESs above the 95th percentile
            q95<-as.numeric(quantile(abs(null_es),0.95))
            fit<-pareto.fit(abs(null_es),threshold=q95)
            newp.value<-ppareto(abs(es), threshold=q95,
                                exponent=fit$exponent, lower.tail=FALSE)/20
            # If Pareto cannot infer small ESs take the permutation p-value
            if(is.na(newp.value)){
                newp.value<-p.value
            }
            p.value<-newp.value
        }
    }

    # Calculate the normalized enrichment score
    nes<-p2z(p.value)*sign(es)

    gsea.obj<-list(
        es=es,
        nes=nes,
        p.value=p.value,
        ledge=ledge_names,
        running_score=running_score,
        set=set,
        reflist=reflist,
        inSet=inSet
    )
    class(gsea.obj)<-"gsea"
    return(gsea.obj)
}


# Calculate Null Distribution for GSEA
null_gsea<-function(set,reflist,w=1,np=1000){
    gsea_null <- rep(0, np)
    gsea_null <- sapply(1:np, function(i) {
        # Identify indexes of set within the sorted reference list
        inSet <- rep(0, length(reflist))
        inSet[which(names(reflist) %in% set)] <- 1

        # By sampling the order of the set elements, we get the real permutation
        null_inSet <- inSet[sample(1:length(inSet))]

        # Same as before, cumulative sums of hits and nonhits
        null_hit<-abs(reflist*null_inSet)
        null_hit<-null_hit^w
        null_hit <- cumsum(null_hit)
        null_hit <- null_hit/null_hit[length(null_hit)]
        null_miss <- cumsum(1-null_inSet)
        null_miss <- null_miss/null_miss[length(null_miss)]
        # And dependending on the cumulative sums, null running sum and
        # null enrichment score
        null_running_score <- null_hit - null_miss

        # The ES is just he maximum or the minimum
        if(abs(max(null_running_score))>abs(min(null_running_score))){
            null_es<-max(null_running_score)
        } else {
            null_es<-min(null_running_score)
        }
        return(null_es)
    })
    class(gsea_null)<-"gsea_null"
    return(gsea_null)
}

#' Plot GSEA results
#'
#' This function generates a GSEA plot from a gsea object
#'
#' @param gsea.obj GSEA object produced by the \code{gsea} function
#' @param twoColors the two colors to use for positive[1] and negative[2]
#' enrichment scores
#' @param colBarcode The color of the barcode
#' @param plotNames Logical. Should the set names be plotted?
#' @param title String to be plotted above the Running Enrichment Score
#' @param bottomYlabel String for the Y label of the bottom plot
#' @param bottomTitle String for the title of the bottom part of the plot
#' @param ext_nes Provide a NES from an external calculation
#' @param omit_middle If TRUE, will not plot the running score
#' (FALSE by default)
#' @return Nothing, a plot is generated in the default output device
#' @examples
#' reflist<-setNames(-sort(rnorm(1000)),paste0('gene',1:1000))
#' set<-paste0('gene',sample(1:200,50))
#' obj<-gsea(reflist,set,method='pareto',np=1000)
#' plot_gsea(obj)
#' @export
plot_gsea <- function(gsea.obj, twoColors = c("red",
                                              "blue"),
                      plotNames = FALSE, colBarcode = "black",
                      title = "Running Enrichment Score",
                      bottomTitle = "List Values",
                      bottomYlabel = "Signature values",
                      ext_nes = NULL, omit_middle = FALSE) {
    # Extract parameters from the gsea object
    es <- gsea.obj$es
    nes <- gsea.obj$nes
    p.value <- gsea.obj$p.value
    ledge <- gsea.obj$ledge
    running_score <- gsea.obj$running_score
    set <- gsea.obj$set
    reflist <- gsea.obj$reflist
    inSet <- gsea.obj$inSet

    # Define plot borders? Who wrote the
    # original code has a non-euclidean mind
    min.RES <- min(running_score)
    max.RES <- max(running_score)
    delta <- (max.RES - min.RES) * 0.5
    min.plot <- min.RES
    max.plot <- max.RES
    max.corr <- max(reflist)
    min.corr <- min(reflist)
    corr0.line <- (-min.corr/(max.corr -
                                  min.corr)) * 1.25 * delta + min.plot

    if (es < 0) {
        l.ledge.ref.plot <- length(reflist) -
            length(ledge)
    } else {
        l.ledge.ref.plot <- length(ledge)
    }

    # Define colors (red is positive nes,
    # blue is negative nes)
    if (nes > 0) {
        col.f <- twoColors[1]
    } else {
        col.f <- twoColors[2]
    }
    N <- length(reflist)
    ind <- 1:N


    ### Define layout, let's end this
    ### destructive putting everything in a
    ### single window
    if (omit_middle) {
        layoutMatrix <- rbind(1, 2)
        layout(layoutMatrix, heights = c(1,
                                         2))
    } else {
        layoutMatrix <- rbind(1, 2, 3)
        layout(layoutMatrix, heights = c(1,
                                         4, 2))
    }


    # layout.show(n=3)

    ### PLOT 1: barcode-like enrichment tags
    par(mar = c(0, 4.1, 2, 2.1))
    plot(0, col = "white", xlim = c(1, N),
         ylim = c(0, 10), xaxt = "n", yaxt = "n",
         type = "n", frame.plot = FALSE, xlab = "",
         ylab = "", xaxs = "r", yaxs = "r",
         main = paste("Number of elements: ",
                      N, " (in full list), ", length(set),
                      " (in element set)", sep = "",
                      collapse = ""))
    for (position in 1:N) {
        if (inSet[position] == 1) {
            if (N < 50 & length(set) <= 10) {
                rect(xleft = position - 0.2,
                     ybottom = 0, xright = position +
                         0.2, ytop = 10, col = colBarcode,
                     border = NA)
            } else {
                abline(v = position, lwd = 2,
                       col = colBarcode)
            }
            if (plotNames) {
                text(labels = names(reflist[position]),
                     x = position - 0.2, y = 0,
                     srt = 90, offset = 0, pos = 4,
                     font = 2)
            }
        }
    }

    ### PLOT 2: The running sum plot
    if (!omit_middle) {
        par(mar = c(2, 4.1, 2, 2.1))
        plot(ind, running_score, sub = "",
             xlab = "", ylab = "Enrichment Score",
             xlim = c(1, N), ylim = c(min.plot,
                                      max.plot), type = "l", pch = 20,
             lwd = 4, cex = 1, col = col.f,
             xaxs = "r", yaxs = "r", main = title)
        grid(col = "dark grey", lty = 2)

        # This is important: zero running score
        # line
        lines(c(1, N), c(0, 0), lwd = 1,
              lty = 1, cex = 1)

        # This is also important: it's the max
        # enrichment vertical line (aka LEADING
        # EDGE)
        lines(c(l.ledge.ref.plot, l.ledge.ref.plot),
              c(min.plot, max.plot), lwd = 1,
              lty = 1, cex = 1)

        if (es >= 0) {
            legend_position <- "topright"
        } else {
            legend_position <- "topleft"
        }

        # If an external NES is not provided, the
        # standard GSEA one is shown
        if (is.null(ext_nes)) {
            legend(legend_position,
                   legend = c(paste("ES = ",
                                    signif(es, 3), sep = ""),
                              paste("NES = ", signif(nes,
                                                     3), sep = ""),
                              paste("p-value = ",
                                    signif(p.value, 3), sep = "")),
                   bg = "white")
        } else {
            legend(legend_position,
                   legend = c(paste("NES = ",
                                    signif(ext_nes, 3), sep = ""),
                              paste("p-value = ",
                                    signif(z2p(ext_nes),
                                           3), sep = "")), bg = "white")
        }
    }


    ### Plot3: weight values
    if (omit_middle) {
        bottomMain <- title
    } else {
        bottomMain <- bottomTitle
    }

    par(mar = c(2, 4.1, 2, 2.1))
    plot(ind, reflist, type = "l", pch = 20,
         lwd = 3, xlim = c(1, N), cex = 1,
         col = 1, xaxs = "r", yaxs = "r",
         main = bottomMain, ylab = bottomYlabel,
         cex.axis = 0.8)
    grid(col = "dark grey", lty = 2)
    # zero correlation horizontal line
    lines(c(1, N), c(corr0.line, corr0.line),
          lwd = 1, lty = 1, cex = 1, col = 1)

    if (omit_middle) {
        # If an external NES is not provided, the
        # standard GSEA one is shown
        if (is.null(ext_nes)) {
            legend("top",
                   legend = c(paste("ES = ",
                                    signif(es, 3), sep = ""),
                              paste("NES = ",
                                    signif(nes,
                                           3), sep = ""),
                              paste("p-value = ",
                                    signif(p.value, 3), sep = "")),
                   bg = "white")
        } else {
            legend("top",
                   legend = c(paste("NES = ",
                                    signif(ext_nes, 3), sep = ""),
                              paste("p-value = ", signif(z2p(ext_nes),
                                                         3), sep = "")),
                   bg = "white")
        }
    }

}



# Probability density of Pareto distributions
dpareto <- function(x, threshold = 1, exponent,
                    log = FALSE) {
    # Avoid doing limited-precision
    # arithmetic followed by logs if we want
    # the log!
    if (!log) {
        prefactor <- (exponent - 1)/threshold
        f <- function(x) {
            prefactor * (x/threshold)^(-exponent)
        }
    } else {
        prefactor.log <- log(exponent - 1) -
            log(threshold)
        f <- function(x) {
            prefactor.log - exponent * (log(x) -
                                            log(threshold))
        }
    }
    d <- ifelse(x < threshold, NA, f(x))
    return(d)
}

# Cumulative distribution function of the Pareto distributions
ppareto <- function(x, threshold = 1, exponent,
                    lower.tail = TRUE) {
    if (!lower.tail) {
        f <- function(x) {
            (x/threshold)^(1 - exponent)
        }
    }
    if (lower.tail) {
        f <- function(x) {
            1 - (x/threshold)^(1 - exponent)
        }
    }
    p <- ifelse(x < threshold, NA, f(x))
    return(p)
}

# Estimate parameters of Pareto distribution
pareto.fit <- function(data, threshold) {
    data <- data[data >= threshold]
    n <- length(data)
    x <- data/threshold
    alpha <- 1 + n/sum(log(x))
    # Calculate Log-Likelihood
    loglike <- sum(dpareto(data, threshold = threshold,
                           exponent = alpha, log = TRUE))
    # KS distance
    newdata <- data[data >= threshold]
    d <- suppressWarnings(ks.test(newdata,
                                  ppareto, threshold = threshold,
                                  exponent = alpha))
    ks.dist <- as.vector(d$statistic)
    fit <- list(type = "pareto", exponent = alpha,
                xmin = threshold, loglike = loglike,
                ks.dist = ks.dist, samples.over.threshold = n)
    return(fit)
}


#' z2p
#'
#' This function gives a gaussian p-value corresponding to the provided Z-score
#'
#' @param z a Z score
#' @return a p-value
#' @examples
#' z<-1.96
#' z2p(z)
#' @export
z2p <- function(z) {
    pnorm(abs(z), lower.tail = FALSE) * 2
}

#' p2z
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param p a p-value
#' @return z a Z score
#' @examples
#' p<-0.05
#' p2z(p)
#' @export
p2z <- function(p) {
    qnorm(p/2, lower.tail = FALSE)
}

#' Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,3,5,2,3)
#' stouffer(zs)
#' @export
stouffer <- function(x) {
    Z <- sum(x)/sqrt(length(x))
    return(Z)
}

#' Weighted Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @param w weight for each Z score
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,-3,5,2,3)
#' ws<-c(1,10,1,2,1)
#' wstouffer(zs,ws)
#' @export
wstouffer <- function(x, w) {
    Z <- sum(x * w)/sqrt(sum(w^2))
    return(Z)
}


#' Fisher integration of p-values
#'
#' This function applies the Fisher integration of pvalues
#'
#' @param ps a vector of p-values
#' @return p.val an integrated p-value
#' @examples
#' ps<-c(0.01,0.05,0.03,0.2)
#' fisherp(ps)
#' @export
fisherp <- function(ps) {
    Xsq <- -2 * sum(log(ps))
    p.val <- pchisq(Xsq, df = 2 * length(ps),
                    lower.tail = FALSE)
    # p<-c(Xsq = Xsq, p.value = p.val)
    return(p.val)
}


#' Slice
#'
#' This function prints a slice of a matrix
#'
#' @param matrix A matrix
#' @return prints it
#' @examples
#' set.seed(1)
#' example<-matrix(rnorm(1000),nrow=100,ncol=10)
#' slice(example)
#' @export
slice <- function(matrix) {
    if (nrow(matrix) < 5) {
        stop("Input matrix has less than 5 rows")
    }
    if (ncol(matrix) < 5) {
        stop("Input matrix has less than 5 columns")
    }
    print(matrix[1:5, 1:5])
}

#' 2-way GSEA
#' GSEA Gene set enrichment analysis of two complementary gene sets using gsea
#' @param reflist named vector of reference scores
#' @param set1 element set 1
#' @param set2 element set 1
#' @param method one of 'permutation' or 'pareto'
#' @param w exponent used to raise the supplied scores. Default is 1 (original
#' scores unchanged)
#' @param np Number of permutations (Default: 1000)
#' @param gsea_null a GSEA null distribution (Optional)
#' @return A list of 2 GSEA objects. Each of which is a list of components:
#' \describe{
#' \item{ES}{The enrichment score}
#' \item{NES}{The normalized enrichment socre}
#' \item{ledge}{The items in the leading edge}
#' \item{p.value}{The permutation-based p-value}
#' }
#' @examples
#' reflist<-setNames(-sort(rnorm(1000)),paste0('gene',1:1000))
#' set1<-paste0('gene',sample(1:200,50))
#' set2<-paste0('gene',sample(801:1000,50))
#' obj<-gsea2(reflist,set1,set2,method='pareto',np=1000)
#' obj$p.value
#' @export
gsea2 <- function (reflist, set1,set2, method = c("permutation", "pareto"),
                   np = 1000, w = 1, gsea_null = NULL) {
    g1 <- gsea(reflist,set1,method=method,np = np, w = w, gsea_null = gsea_null)
    g2 <- gsea(reflist,set2,method=method,np = np, w = w, gsea_null = gsea_null)
    ix <- order(reflist, decreasing = T)
    reflist <- reflist[ix]

    gsea.obj <- list(
        es1 = g1$es,es2 = g2$es,
        nes1 = g1$nes, nes2 =  g2$nes,
        p.value1 = g1$p.value,p.value2 =  g2$p.value,
        ledge1 = g1$ledge,ledge2 =  g2$ledge,
        running_score1 = g1$running_score,running_score2 =  g2$running_score,
        set1 = set1,set2 = set2,
        reflist = reflist,
        inSet1 = g1$inSet,inSet2 = g2$inSet)
    return(gsea.obj)
}

#' Plot 2-way GSEA results
#'
#' This function generates a GSEA plot from a gsea object
#'
#' @param gsea.obj GSEA object produced by the \code{gsea} function
#' @param twoColors the two colors to use for positive[1] and negative[2]
#' enrichment scores, and of the barcodes
#' @param plotNames Logical. Should the set names be plotted?
#' @param title String to be plotted above the Running Enrichment Score
#' @param bottomTitle String for the title of the bottom part of the plot
#' @param bottomYlabel String for the Y label of the bottom plot
#' (FALSE by default)
#' @return Nothing, a plot is generated in the default output device
#' @examples
#' reflist<-setNames(-sort(rnorm(1000)),paste0('gene',1:1000))
#' set1<-paste0('gene',sample(1:200,50))
#' set2<-paste0('gene',sample(801:1000,50))
#' obj<-gsea2(reflist,set1,set2,method='pareto',np=1000)
#' plot_gsea2(obj)
#' @export
plot_gsea2 <- function (gsea.obj, twoColors = c("red", "blue"),
                        plotNames = FALSE,
                        title = "Running Enrichment Score",
                        bottomTitle = "List Values",
                        bottomYlabel = "Signature values") {
    es1 <- gsea.obj$es1; es2 <- gsea.obj$es2;
    nes1 <- gsea.obj$nes1;nes2 <- gsea.obj$nes2;
    p.value1 <- gsea.obj$p.value1; p.value2 <- gsea.obj$p.value2;
    ledge1 <- gsea.obj$ledge1; ledge2 <- gsea.obj$ledge2;
    running_score1 <- gsea.obj$running_score1;running_score2 <- gsea.obj$running_score2;
    set1 <- gsea.obj$set1;set2 <- gsea.obj$set2;
    inSet1 <- gsea.obj$inSet1; inSet2 <- gsea.obj$inSet2;
    reflist <- gsea.obj$reflist

    min.plot <- min(running_score1,running_score2)
    max.plot <- max(running_score1,running_score2)

    if (es1 < 0) {
        l.ledge.ref.plot1 <- length(reflist) - length(ledge1)
    }else {
        l.ledge.ref.plot1 <- length(ledge1)
    }

    if (es2 < 0) {
        l.ledge.ref.plot2 <- length(reflist) - length(ledge2)
    }else {
        l.ledge.ref.plot2 <- length(ledge2)
    }

    N <- length(reflist)
    ind <- 1:N

    layoutMatrix <- rbind(1, 2, 3, 4)
    layout(layoutMatrix, heights = c(1, 5, 1, 2))
    ##-----------------------------------------
    par(mar = c(0, 4.1, 2, 2.1))
    plot(0, col = "white", xlim = c(1, N), ylim = c(0, 10), xaxt = "n",
         yaxt = "n", type = "n", frame.plot = FALSE, xlab = "",
         ylab = "", xaxs = "r", yaxs = "r", main = paste("Number of elements: ",
                                                         N, " (in full list), ", length(set1), " (in element set1)",
                                                         sep = "", collapse = ""))
    for (position in 1:N) {
        if (inSet1[position] == 1) {
            if (N < 50 | length(set1) <= 10) {
                rect(xleft = position - 0.2, ybottom = 0, xright = position +
                         0.2, ytop = 10, col = twoColors[1], border = NA)
            }
            else {
                abline(v = position, lwd = 1, col = twoColors[1])
            }
            if (plotNames) {
                text(labels = names(reflist[position]), x = position -
                         0.2, y = 0, srt = 90, offset = 0, pos = 4,
                     font = 2)
            }
        }
    }
    ##--------------------------------
    par(mar = c(0, 4.1, 2, 2.1))
    plot(x=NULL,y=NULL,xlim=range(c(1,N)), ylim=range(c(min.plot, max.plot)), sub = "", xlab = "", ylab = "Running Enrichment Score",
         pch = 20, lwd = 2, cex = 1, xaxt = "n",
         yaxs = "r", main = title)
    grid(col = "light grey", lty = 2)
    lines(ind, running_score1, sub = "", xlab = "", ylab = "Enrichment Score",
          lwd = 2, cex = 1, col = twoColors[1])
    lines(ind, running_score2, sub = "", xlab = "", ylab = "Enrichment Score",
          lwd = 2, cex = 1, col = twoColors[2])
    lines(c(1, N), c(0, 0), lwd = 1, lty = 1, cex = 1)
    lines(c(l.ledge.ref.plot1, l.ledge.ref.plot1), c(min(running_score1),
                                                     max(running_score1)), lwd = 1, lty = 1, cex = 1)
    lines(c(l.ledge.ref.plot2, l.ledge.ref.plot2), c(min(running_score2),
                                                     max(running_score2)), lwd = 1, lty = 1, cex = 1)
    if(es1 >=0){legside1="topright"
    }else{legside1="topleft"}
    legend(legside1, legend = c(paste("NES1 = ", signif(nes1, 3), sep = ""),
                                paste("P1 = ", sprintf("%.2e",p.value1), sep = "")),
           bg = "white",box.col=twoColors[1],cex=1.5)
    if(es2 >=0){legside2="bottomright"
    }else{legside2="bottomleft"}
    legend(legside2, legend = c(paste("NES2 = ", signif(nes2, 3), sep = ""),
                                paste("P2 = ", sprintf("%.2e",p.value2), sep = "")),
           bg = "white",box.col=twoColors[2],cex=1.5)
    # Integration
    p.value.int<-fisherp(c(p.value1,p.value2))
    legend("center",legend=paste0("pFisher = ",sprintf("%.2e",p.value.int)),cex=1.5)

    ##----------------------------------
    par(mar = c(0, 4.1, 2, 2.1))
    plot(0, col = "white", xlim = c(1, N), ylim = c(0, 10), xaxt = "n",
         yaxt = "n", type = "n", frame.plot = FALSE, xlab = "",
         ylab = "", xaxs = "r", yaxs = "r", main = paste("Number of elements: ",
                                                         N, " (in full list), ", length(set2), " (in element set2)",
                                                         sep = "", collapse = ""))
    for (position in 1:N) {
        if (inSet2[position] == 1) {
            if (N < 50 | length(set2) <= 10) {
                rect(xleft = position - 0.2, ybottom = 0, xright = position +
                         0.2, ytop = 10, col = twoColors[2], border = NA)
            }
            else {
                abline(v = position, lwd = 1, col = twoColors[2])
            }
            if (plotNames) {
                text(labels = names(reflist[position]), x = position -
                         0.2, y = 0, srt = 90, offset = 0, pos = 4,
                     font = 2)
            }
        }
    }

    #-----------------------------------
    par(mar = c(2, 4.1, 2, 2.1))
    plot(x=NULL,y=NULL,xlim=range(c(0,N)),
         ylim=range(c(min(reflist),max(reflist))), type = "b", lwd = 2,
         cex = 1, col = 1, xaxs = "r", yaxs = "r",main=bottomTitle,
         ylab=bottomYlabel)
    grid(col = "light grey", lty = 2)
    abline(h=0, lwd = 1,lty = 1, cex = 1, col = 1)
    lines(ind,reflist,lty=1,lwd=2)
}


