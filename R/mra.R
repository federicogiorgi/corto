#' Perform Master Regulator Analysis (mra).
#'
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
#' Default is 1000 if expmat2 is provided, and 10 if expmat2 is not provided (single sample mra).
#' @param nthreads The number of threads to use for generating null signatures. Default is 1
#' @param atacseq An optional 3 column matrix derived from an ATAC-Seq analysis, indicating
#' 1) gene symbol, 2) -log10(FDR)*sing(log2FC) of an ATAC-Seq design, 3) distance from TSS.
#' If provided, the output will contain an _atacseq_ field.
#' @param verbose Boolean, whether to print full messages on progress analysis. Default is FALSE
#' @return A list summarizing the master regulator analysis
#' \itemize{
#' \item nes: the normalized enrichment score: positive if the centroid/TF network is upregulated
#' in expmat1 vs expmat2 (or in expmat1 vs the mean of the dataset), negative if downregulated. A
#' vector in multisample mode, a matrix in sample-by-sample mode.
#' \item pvalue: the pvalue of the enrichment.
#' \item sig: the calculated signature (useful for plotting).
#' \item regulon: the original regulon used in the analysis (but filtered for _minsize_)
#' \item atac: Optionally present if atacseq data is provided. For each centroid/TF a number
#' ranging from 0 to 1 will indicate the fraction of changes in activity due to promoter effects
#' rather than distal effects.
#' }
#' @export
mra<-function(expmat1,expmat2=NULL,regulon,minsize=10,nperm=NULL,nthreads=2,verbose=FALSE,
              atacseq=NULL){
    # Setting default nperm if not set by the user
    if(is.null(nperm)){
        if(is.null(expmat2)){
            nperm<-10
        } else{
            nperm<-1000
        }
    }

    # Filter by minsize
    regsizes<-sapply(regulon,function(x){length(x$likelihood)})
    regulon<-regulon[regsizes>=minsize]
    regsizes<-regsizes[names(regulon)]

    # First we create a centroid/target matrix
    centroids<-sort(names(regulon))
    targets<-sort(unique(unlist(sapply(regulon,function(x){names(x$tfmode)}))))
    netmat<-matrix(0,nrow=length(centroids),ncol=length(targets),dimnames=list(centroids,targets))
    # Then we fill it
    for(centroid in centroids){
        vec<-sign(regulon[[centroid]]$tfmode)*regulon[[centroid]]$likelihood
        netmat[centroid,names(vec)]<-vec
    }

    # Case 1: expmat2 is provided as control
    if(!is.null(expmat2)){
        # We create a signature
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

        # # Then we normalize by regulon length, otherwise TFs with small networks will be penalized
        # netmat<-t(apply(netmat,1,function(x){10*x/sum(abs(x))}))
        # netmat[is.na(netmat)]<-0

        # Calculate unnormalized scores
        scores<-(netmat%*%sig)[,1]

        # # Check Correlation of scores with regsizes
        # regsizes<-regsizes[names(scores)]
        # plot(regsizes,scores)
        # mtext(paste0("Cor=",cor(regsizes,scores)))
        # # Interestingly here there is no correlation

        #cands<-c("BAZ1A","HCFC1","HDGF","MAZ","ZNF146","ZNF532")
        #scores[cands]

        # Permuted signatures
        nullsigperm1<-function(seed=0,expmat1,expmat2,netmat){
            permat<-cbind(expmat1,expmat2)
            set.seed(seed)
            permat<-permat[sample(nrow(permat)),]
            permat<-permat[,sample(ncol(permat))]
            nullsig<-setNames(apply(permat,1,function(x){
                if(!is.null(expmat2)){
                    x1<-x[1:ncol(expmat1)]
                    x2<-x[(ncol(expmat1)+1):length(x)]
                    tt<-t.test(x1,x2)
                } else {
                    x<-x-mean(x)
                    tt<-t.test(x)
                }
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
                             FUN=nullsigperm1,
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
            if(p==0){p<-.Machine$double.xmin}
            mynes<-p2z(p)*sign(myscore)
            if(myscore==0){mynes<-0}
            return(mynes)
        })
        outlist<-list(nes=nes,pvalue=z2p(nes),sig=sig,regulon=regulon)
        return(outlist)

        # # Check correlation regsize nes
        # regsizes<-regsizes[names(nes)]
        # plot(regsizes,nes,col="grey")
        # mtext(paste0("Cor=",cor(regsizes,nes)))
        # text(regsizes[cands],nes[cands],labels=cands,col="black",cex=2,font=2)
        #
        # # Check correlation scores nes
        # scores<-scores[names(nes)]
        # plot(scores,nes)
        # mtext(paste0("Cor=",cor(scores,nes)))

    } else { # Case 2: expmat2 is not provided
        # Create signatures
        sigmat<-t(apply(expmat1,1,function(x){
            y<-(x-mean(x))/sd(x)
            return(y)
        }))
        # sigs<-t(scale(t(expmat1))) # Exactly the same
        sigmat[is.na(sigmat)]<-0
        common<-intersect(rownames(sigmat),colnames(netmat))
        sigmat<-sigmat[common,]
        netmat<-netmat[,common]

        # Then we slightly reduce the contribution of targets shared by multiple centroids
        # But only targets in the signature
        netmat<-apply(netmat,2,function(x){x/(sum(x!=0)^0.5)})

        # Calculate unnormalized scores
        scores<-(netmat%*%sigmat)

        # Permuted signatures
        nullsigperm2<-function(seed=0,expmat1,netmat){
            permat<-expmat1
            set.seed(seed)
            permat<-permat[sample(nrow(permat)),]
            permat<-permat[,sample(ncol(permat))]
            nullsigmat<-t(apply(permat,1,function(x){
                y<-(x-mean(x))/sd(x)
                return(y)
            }))
            nullsigmat<-nullsigmat[colnames(netmat),]
            nullscores<-(netmat%*%nullsigmat)
            return(nullscores)
        }

        # The pbapply snippet
        cl<-parallel::makeCluster(nthreads)
        #invisible(parallel::clusterExport(cl=cl,varlist=c("nthreads")))
        #invisible(parallel::clusterEvalQ(cl= cl,library(dplyr)))
        nullscores<-pblapply(cl=cl,
                             X=1:nperm,
                             FUN=nullsigperm2,
                             expmat1=expmat1,netmat=netmat
        )
        parallel::stopCluster(cl)
        nullscores<-do.call(cbind,nullscores)

        # Compare scores with nullscores
        # Fit gaussians
        nes<-t(apply(cbind(scores,nullscores),1,function(x){
            myscores<-x[1:ncol(scores)]
            morescores<-x[(ncol(scores)+1):length(x)]
            mu<-mean(morescores)
            sigma<-sd(morescores)
            ps<-pnorm(abs(myscores),mean=mu,sd=sigma,lower.tail=FALSE)*2
            ps[ps==0]<-.Machine$double.xmin
            mynes<-p2z(ps)*sign(myscores)
            mynes[myscores==0]<-0
            return(mynes)
        }))
        return(nes)
        # pnb<-nes["MYCN",names(hashtags)[hashtags=="pnb"]]
        # ctr<-nes["MYCN",names(hashtags)[hashtags=="ctr"]]
        # plot(density(pnb),col="red",lwd=3,main="MYCN Activity")
        # lines(density(ctr),col="blue",lwd=3)
        # legend("topright",col=c("red","blue"),legend=c("Panobinostat","Control"),lwd=3)
        # abline(v=0,lty=2)
    }
}

#' Plot a master regulator analysis
#'
#' Plotting function for master regulator analysis performed by the _mra_ function
#' @param mraobj The input object, output of the function mra
#' @param mrs Either a numeric value indicating how many MRs to show, sorted by
#' significance, or a character vector specifying which TFs to show
#' @return A plot is generated
#' @export
mraplot<-function(mraobj,mrs=NULL){
    # Checks ----
    if(is.numeric(mrs)){
        mrs<-names(sort(abs(mraobj$nes),decreasing=TRUE))[1:mrs]
    } else if(is.character(mrs)){
        mrs<-mrs
    } else {
        stop("Provide mrs as a number to show top mrs, or a vector of mrs")
    }

    # Define layout ----
    ltitle<-c(1,1,1,2,2)
    lblocks<-matrix(NA,nrow=0,ncol=length(ltitle))
    for(i in 1:length(mrs)){
        base<-3+(6*(i-1))
        lblock<-rbind(
            c(base,base+1,base+2,base+3,base+3),
            c(base+4,base+4,base+4,base+3,base+3),
            c(base+5,base+5,base+5,base+3,base+3)
        )
        lblocks<-rbind(lblocks,lblock)
    }
    lmatrix<-rbind(
        ltitle,
        lblocks
    )

    # Start plotting ----
    layout(lmatrix,heights=c(1,rep(2,nrow(lmatrix)-1)))
    # layout.show(n=max(lmatrix))

    # Function for plotting text only ----
    titplot<-function(title,box=TRUE,cex=2,bgcol="white",bold=FALSE){
        if(bold){font<-2}else{font<-1}
        opar<-par()$mar
        par(mar = c(0.1,0.1,0.1,0.1))
        plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
        rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=bgcol)
        text(x=0.5,y=0.5,title, cex = cex, col = "black",font=font)
        # if(bottom){
        #     axis(1,lwd=10,tick=FALSE,labels=FALSE)
        # }
        if(box){
            box()
        }
        par(mar=opar)
    }
    # Panels 1 and 2 are titles ----
    titplot("corto - Master Regulator Analysis")
    titplot("Top targets")

    # Prefetch signature ----
    ranksig<-rank(mraobj$sig)
    ranksiginv<-rank(-mraobj$sig)
    ranksigabs<-rank(abs(mraobj$sig))


    # Fill the plot with MR blocks ----
    for(mr in mrs){
        # Name of the MR
        titplot(mr,cex=4)

        ### NES ----
        bgcol<-"white"
        if(mraobj$nes[mr]>0&mraobj$pvalue[mr]<=0.01){
            bgcol<-"salmon"
        } else if(mraobj$nes[mr]<0&mraobj$pvalue[mr]<=0.01){
            bgcol<-"cornflowerblue"
        }
        titplot(paste0("NES=",round(mraobj$nes[mr],2)),bgcol=bgcol)
        ### p-value ----
        bold<-FALSE
        if(mraobj$pvalue[mr]<=0.01){bold<-TRUE}
        titplot(paste0("p=",signif(mraobj$pvalue[mr],3)),bold=bold)

        ### Network ----
        # We show maybe the top 12 targets
        # titplot("Network goes here")
        opar<-par()$mar
        par(mar=c(0,0,0,0),xaxs="i")
        plot(0,xlim=c(-1.5,1.5),ylim=c(-1.2,1.2),type="n",ann=F,bty='n',xaxt='n',yaxt='n')

        # Which targets to show
        targets<-names(mraobj$regulon[[mr]]$tfmode)
        targets<-intersect(targets,names(ranksig))
        toshow<-names(sort(ranksigabs[targets],decreasing=TRUE)[1:12])

        # Colors
        col<-rep("white",12)
        col[mraobj$sig[toshow]<0]<-"#6495ED66" # cornflowerblue
        col[mraobj$sig[toshow]>0]<-"#FF8C6966" # salmon

        # Type of regulation (mode of action)
        moa<-rep("unknown",12)
        tfmodehere<-mraobj$regulon[[mr]]$tfmode[toshow]
        moa[tfmodehere>0]<-"activation"
        moa[tfmodehere<0]<-"repression"
        angle<-moa
        angle[moa=="activation"]<-30
        angle[moa=="repression"]<-90
        angle[moa=="unknown"]<-0

        # Circles Clock
        draw.circle(0,1,0.1,border="gray90",col=col[1])
        arrows(0,0.1,0,0.9,lwd=2,col="gray50",angle=angle[1])
        text(0,1,font=2,cex=1.5,labels=toshow[1])

        draw.circle(sin(pi/6),cos(pi/6),0.1,border="gray90",col=col[2])
        arrows(0,0.1,sin(pi/6)-0.1,cos(pi/6)-0.1,lwd=2,col="gray50",angle=angle[2])
        text(sin(pi/6),cos(pi/6),font=2,cex=1.5,labels=toshow[2])

        draw.circle(sin(pi/3),cos(pi/3),0.1,border="gray90",col=col[3])
        arrows(0,0.1,sin(pi/3)-0.2,cos(pi/3)-0.1,lwd=2,col="gray50",angle=angle[3])
        text(sin(pi/3),cos(pi/3),font=2,cex=1.5,labels=toshow[3])

        draw.circle(1,0,0.1,border="gray90",col=col[4])
        arrows(0.1,0,0.7,0,lwd=2,col="gray50",angle=angle[4])
        text(1,0,font=2,cex=1.5,labels=toshow[4])

        draw.circle(sin(pi/3),-cos(pi/3),0.1,border="gray90",col=col[5])
        arrows(0,-0.1,sin(pi/3)-0.2,-cos(pi/3)+0.1,lwd=2,col="gray50",angle=angle[5])
        text(sin(pi/3),-cos(pi/3),font=2,cex=1.5,labels=toshow[5])

        draw.circle(sin(pi/6),-cos(pi/6),0.1,border="gray90",col=col[6])
        arrows(0,-0.1,sin(pi/6)-0.1,-cos(pi/6)+0.1,lwd=2,col="gray50",angle=angle[6])
        text(sin(pi/6),-cos(pi/6),font=2,cex=1.5,labels=toshow[6])

        draw.circle(0,-1,0.1,border="gray90",col=col[7])
        arrows(0,-0.1,0,-0.9,lwd=2,col="gray50",angle=angle[7])
        text(0,-1,font=2,cex=1.5,labels=toshow[7])

        draw.circle(-sin(pi/6),-cos(pi/6),0.1,border="gray90",col=col[8])
        arrows(0,-0.1,-sin(pi/6)+0.1,-cos(pi/6)+0.1,lwd=2,col="gray50",angle=angle[8])
        text(-sin(pi/6),-cos(pi/6),font=2,cex=1.5,labels=toshow[8])

        draw.circle(-sin(pi/3),-cos(pi/3),0.1,border="gray90",col=col[9])
        arrows(0,-0.1,-sin(pi/3)+0.2,-cos(pi/3)+0.1,lwd=2,col="gray50",angle=angle[9])
        text(-sin(pi/3),-cos(pi/3),font=2,cex=1.5,labels=toshow[9])

        draw.circle(-1,0,0.1,border="gray90",col=col[10])
        arrows(-0.1,0,-0.7,0,lwd=2,col="gray50",angle=angle[10])
        text(-1,0,font=2,cex=1.5,labels=toshow[10])

        draw.circle(-sin(pi/3),cos(pi/3),0.1,border="gray90",col=col[11])
        arrows(0,0.1,-sin(pi/3)+0.1,cos(pi/3)-0.1,lwd=2,col="gray50",angle=angle[11])
        text(-sin(pi/3),cos(pi/3),font=2,cex=1.5,labels=toshow[11])

        draw.circle(-sin(pi/6),cos(pi/6),0.1,border="gray90",col=col[12])
        arrows(0,0.1,-sin(pi/6)+0.2,cos(pi/6)-0.1,lwd=2,col="gray50",angle=angle[12])
        text(-sin(pi/6),cos(pi/6),font=2,cex=1.5,labels=toshow[12])


        # Plot Centroid
        draw.circle(0,0,0.1,border="gray90",col="#00000011")
        text(0,0,labels=mr,cex=2,font=2)

        box()
        par(mar=opar)



        ### Signature ----
        # Define Transparency for barcode plot coloring
        if(mraobj$nes[mr]>0&mraobj$pvalue[mr]<=0.01){
            transp<-255*(ranksig^3)/(length(ranksig)^3)
            transp<-as.hexmode(round(transp))
            transp<-format(transp,width=2)
            names(transp)<-names(ranksig)
            transpinv<-255*(ranksiginv^3)/(length(ranksiginv)^3)
            transpinv<-as.hexmode(round(transpinv))
            transpinv<-format(transpinv,width=2)
            names(transpinv)<-names(ranksig)
        } else if(mraobj$nes[mr]<0&mraobj$pvalue[mr]<=0.01){
            transp<-255*(ranksiginv^3)/(length(ranksiginv)^3)
            transp<-as.hexmode(round(transp))
            transp<-format(transp,width=2)
            names(transp)<-names(ranksig)
            transpinv<-255*(ranksig^3)/(length(ranksig)^3)
            transpinv<-as.hexmode(round(transpinv))
            transpinv<-format(transpinv,width=2)
            names(transpinv)<-names(ranksig)
        } else {
            transp<-255*(ranksigabs^3)/(length(ranksigabs)^3)
            transp<-as.hexmode(round(transp))
            transp<-format(transp,width=2)
            names(transp)<-names(ranksig)
            transpiv<-transp
        }



        # Start plotting
        opar<-par()$mar
        par(mar=c(0,0,0,0),xaxs="i")
        plot(0,xlim=c(1,length(ranksig)),ylim=c(0,1),type="n",xaxt="n")

        targets<-mraobj$regulon[[mr]]$tfmode
        # Positive targets
        postargets<-names(targets[targets>0])
        postargets<-intersect(names(ranksig),postargets)
        pos<-ranksig[postargets]
        if(length(pos)>0){
            segments(pos,0.4,pos,1,col=paste0("#CD0000",transp[postargets]),lwd=3)
        }
        # Negative targets
        negtargets<-names(targets[targets<0])
        negtargets<-intersect(names(ranksig),negtargets)
        neg<-ranksig[negtargets]
        if(length(neg)>0){
            segments(neg,0,neg,0.6,col=paste0("#000080",transpinv[negtargets]),lwd=3)
        }

        # text(0,0.5,labels="-",pos=4,cex=12,font=1)
        # text(length(ranksig),0.5,labels="+",pos=2,cex=10,font=1)
        uparrow<-intToUtf8(8593)
        dnarrow<-intToUtf8(8595)
        text(length(ranksig)/2,0.5,labels=paste0(dnarrow,"   ",uparrow),pos=NULL,cex=5,font=2)
        par(mar=opar)

        ### Keep blank to separate from the next MR
        titplot("")
    }
}

