#' Perform Master Regulator Analysis (mra).
#' The analysis is performed between two groups of samples in the form of expression matrices,
#' with genes/features as rows and samples as columns.
#' @param expmat1 A numeric expression matrix, with genes/features as rows and samples as columns.
#' If only expmat1 is provided (without expmat2), the function will perform a sample-by-sample
#' master regulator analysis, with the mean of the dataset as a reference. If expmat2 is provided,
#' expmat1 will be considered the "treatment" sample set.
#' @param expmat2 A numeric expression matrix, with genes/features as rows and samples as columns.
#' If provided, it will be considered as the "control" or "reference" sample set for expmat1.
#' @param regulon A _regulon_ object, output of the _corto_ function.
#' @param minsize A minimum network size for each centroid/TF to be analyzed. Default is 10.
#' @param nperm The number of times the input data will be permuted to generate null signatures.
#' Default is 100.
#' @param nthreads The number of threads to use for generating null signatures. Default is 1
#' @param atacseq An optional 3 column matrix derived from an ATAC-Seq analysis, indicating
#' 1) gene symbol, 2) -log10(FDR)*sing(log2FC) of an ATAC-Seq design, 3) distance from TSS.
#' If provided, the output will contain an _atacseq_ field.
#' @param verbose Boolean, whether to print full messages on progress analysis
#' @return A list summarizing the master regulator analysis
#' \itemize{
#' \item nes: the normalized enrichment score: positive if the centroid/TF network is upregulated
#' in expmat1 vs expmat2 (or in expmat1 vs the mean of the dataset), negative if downregulated. A
#' vector in multisample mode, a matrix in sample-by-sample mode.
#' \item pvalue: the pvalue of the enrichment.
#' \item sig: the calculated signature (useful for plotting).
#' \item atac: Optionally present if atacseq data is provided. For each centroid/TF a number
#' ranging from 0 to 1 will indicate the fraction of changes in activity due to promoter effects
#' rather than distal effects.
#' }
#' @export
mra<-function(expmat1,expmat2=NULL,regulon,minsize=10,nperm=100,nthreads=2,verbose=FALSE,atacseq=NULL){
    # Filter by minsize TO DO
    regsizes<-sapply(regulon,function(x){length(x$likelihood)})
    regulon<-regulon[regsizes>=minsize]

    # First we create a centroid/target matrix
    centroids<-sort(names(regulon))
    targets<-sort(unique(unlist(sapply(regulon,function(x){names(x$tfmode)}))))
    netmat<-matrix(0,nrow=length(centroids),ncol=length(targets),dimnames=list(centroids,targets))
    # Then we fill it
    for(centroid in centroids){
        vec<-sign(regulon[[centroid]]$tfmode)*regulon[[centroid]]$likelihood
        netmat[centroid,names(vec)]<-vec
    }

    # Multisample: we create a signature
    sig<-setNames(apply(cbind(expmat1,expmat2),1,function(x){
        x1<-x[1:ncol(expmat1)]
        x2<-x[(ncol(expmat1)+1):length(x)]
        tt<-t.test(x1,x2)
        return(tt$statistic)
    }),rownames(expmat1))
    sig<-sig[!is.na(sig)]
    common<-intersect(names(sig),colnames(netmat))
    sig<-sig[common]
    netmat<-netmat[,common]

    # Then we slightly reduce the contribution of targets shared by multiple centroids
    # But only targets in the signature
    netmat2<-apply(netmat,2,function(x){x/(sum(x!=0)^0.5)})
    # plot(netmat[,1:20],netmat2[,1:20])
    # abline(h=0,v=0)
    netmat<-netmat2
    rm(netmat2)

    # Calculate unnormalized scores
    scores<-(netmat%*%sig)[,1]

    #cands<-c("BAZ1A","HCFC1","HDGF","MAZ","ZNF146","ZNF532")
    #scores[cands]

    # Permuted signatures
    nullsigperm<-function(seed=0,expmat1,expmat2,netmat){
        permat<-cbind(expmat1,expmat2)
        set.seed(seed)
        permat<-permat[sample(nrow(permat)),]
        permat<-permat[,sample(ncol(permat))]
        nullsig<-setNames(apply(permat,1,function(x){
            x1<-x[1:ncol(expmat1)]
            x2<-x[(ncol(expmat1)+1):length(x)]
            tt<-t.test(x1,x2)
            return(tt$statistic)
        }),rownames(permat))
        nullsig<-nullsig[colnames(netmat)]
        nullscores<-(netmat%*%nullsig)[,1]
        return(nullscores)
    }

    # The pbapply snippet
    cl<-parallel::makeCluster(nthreads)
    #invisible(parallel::clusterExport(cl=cl,varlist=c("nthreads")))
    #invisible(parallel::clusterEvalQ(cl= cl,library(dplyr)))
    nullscores<-pbsapply(cl=cl,
                         X=1:nperm,
                         FUN=nullsigperm,
                         expmat1=expmat1,expmat2=expmat2,netmat=netmat
    )
    parallel::stopCluster(cl)


    # Compare scores with nullscores
    # Fit gaussians
    nes<-apply(cbind(scores,nullscores),1,function(x){
        myscore<-x[1]
        morescores<-x[2:length(x)]
        mu<-mean(morescores)
        sigma<-sd(morescores)
        p<-pnorm(abs(myscore),mean=mu,sd=sigma,lower.tail=FALSE)*2
        mynes<-p2z(p)*sign(myscore)
        if(myscore==0){mynes<-0}
        return(mynes)
    })

    outlist<-list(nes=nes,pvalue=z2p(nes),sig=sig)
    return(outlist)
}
