library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Matrix)
library(MASS)
#library(biovizBase)


## This is buggy -- does not rearange multiple exons on '-' strand
#transcripts.bad <- extractTranscriptSeqs(yg,cds)
## And this one is extractng exons, not CDS?
##transcripts.bad2 <- extractTranscriptSeqs(yg,txdb)
## Correct, but slow fix

properExtractTranscripts <- function(genome, cds) {
    as(lapply(cds, function(g) if(unique(strand(g))=='+')
                                   unlist(getSeq(genome, g))
                               else
                                   unlist(rev(getSeq(genome, g)))),
       "DNAStringSet")
}

polyaTranscripts <- function(transcripts, n=1) {
    endoapply(transcripts, c, DNAString(paste(rep("AAA", n), collapse="")))
}

## covData <- riboSeqFromBAM(bams, genomeName="sacSer3", txdb=txdb)

## Shifting properly

readPsites <- function(bam, lengths, shifts) {
    stopifnot(length(lengths)==length(shifts))
    aln <- readGAlignments(bams)
    psites <- GRanges()
    for (i in seq_along(length)) {
        l <- lengths[[i]]
        pl <- granges(qnarrow(aln[width(aln)==l], start=l-shifts[[i]], width=1))
        pl$score <- l
        psites <- c(psites, pl)
    }
    return(psites)
}


sitesToWigVector <- function(psites) {
    l <- lapply(seqnames(seqinfo(psites)),
                function(sn) {
                    ps <- psites[seqnames(psites)==sn & strand(psites)=='+']
                    ms <- psites[seqnames(psites)==sn & strand(psites)=='-']
                    tp <- tabulate(start(ps), seqlengths(seqinfo(psites))[sn])
                    tm <- tabulate(start(ms), seqlengths(seqinfo(psites))[sn])
                    list('+'=as(tp, "sparseVector"),
                         '-'=as(tm, "sparseVector"))
                })
    names(l) <- seqnames(seqinfo(psites))
    return(l)
}

sitesToWig <- function(psites) {
    l <- lapply(seqnames(seqinfo(psites)),
                function(sn) {
                    ps <- psites[seqnames(psites)==sn & strand(psites)=='+']
                    ms <- psites[seqnames(psites)==sn & strand(psites)=='-']
                    tp <- tabulate(start(ps))
                    tm <- tabulate(start(ms))
                    c(GRanges(seqnames=sn,
                              IRanges(start=which(tp>0), width=0),
                              strand='+',
                              score=tp[tp>0]),
                      GRanges(seqnames=sn,
                              IRanges(start=which(tm>0), width=0),
                              strand='-',
                              score=tm[tm>0]))
                })
    as(do.call(`c`, l), "UCSCData")
}

exportWigPair <- function(wig, base, name="wig") {
    gru <- as(wig[strand(wig)=='+'], "UCSCData")
    gru@trackLine@name <- paste0(name,'+')
    export.wig(gru, paste0(base,'.plus.wig'))
    gru <- as(wig[strand(wig)=='-'], "UCSCData")
    gru@trackLine@name <- paste0(name,'-')
    export.wig(gru, paste0(base,'.plus.wig'))
}





extractGeneWig <- function(gr, wig, before=0, after=0) {
    ##    chr <- as.character(seqnames(gr))[[1]]  # Assuming all are the same here!
    ##    strand <- as.character(strand(gr))[[1]]
    ## Faster version?
    chr <- as.character(gr@seqnames@values[[1]])  # Assuming all are the same here!
    strand <- as.character(gr@strand@values[[1]])
    w <- wig[[chr]][[strand]]
#    i <<- i+1
#    if (i%%1000==0) print(i)
    if (length(gr)>1)
        pos <- unlist(mapply(seq, start(gr), end(gr)))
    else
        pos <- start(gr):end(gr)
    if (before>0 && strand=='+')
        pos <- c(seq(start(gr)[[1]]-before, start(gr)[[1]]-1), pos)
    if (before>0 && strand=='-')
        pos <- c(pos, seq(tail(end(gr),1)+1, tail(end(gr),1)+before))
    if (after>0 && strand=='+')
        pos <- c(pos, seq(tail(end(gr),1)+1, tail(end(gr),1)+after))
    if (after>0 && strand=='-')
        pos <- c(seq(start(gr)[[1]]-after, start(gr)[[1]]-1), pos)
    res <- as.vector(w[pos])
    if (strand=='-')
        return(rev(res))
    else
        return(res)
}

extractGeneListWigs <- function(grl, wig,...) {
    i <<- 1
    lapply(grl, extractGeneWig, wig,...)
}

extractGeneListWigsbpl <- function(grl, wig,...) {
    require(BiocParallel)
    i <<- 1
    bplapply(grl, extractGeneWig, wig,...)
}

sumGeneWigReadFrame <- function(gw) {
    c(sum(gw[c(TRUE,FALSE,FALSE)]),
      sum(gw[c(FALSE,TRUE,FALSE)]),
      sum(gw[c(FALSE,FALSE,TRUE)]))
}

sumGeneListWigReadFrame <- function(gwl) {
    res <- c(0,0,0)
    for (x in lapply(gwl, sumGeneWigReadFrame))
        res <- res+x
    res
}




doReadFramePlot <- function(fsum,title="") {
    fsum <- fsum/sum(fsum)*1e6
    barplot(fsum, names.arg=c(1,2,3), col=c('red', 'green', 'cyan'),
            ylab='Normalized counts', xlab='Read Frame position',
            main=title, ylim=c(0,6e5))
    text(c(1, 2, 3)*1.2-0.6, y=fsum/2, labels=round(fsum/sum(fsum)*100, 0),
         adj=c(0,0.5))
}



geneCodons <- function(g, before=0, after=0, skipstart=0, skipend=0,win=1) {
    a <- before+skipstart+1
    b <- length(g)-after-skipend
    colSums(matrix(g[a:b], nrow=3*win))
}

codonWigs <- function(gwl) {
    lapply(gwl, geneCodons)
}



geneFitProbability <- function(gw) {
    require(MASS)
    fitdistr(gw, "negative binomial")
}

geneNonPauseLogP <- function(gw, avgcut=0.1, maxcut=4) {
    require(MASS)
    if (mean(gw)>avgcut & max(gw)>maxcut) {
        ## Optimizer tends to check for logP with dnbinom with size=0,
        ##  leading to NANs.  How to deal with it properly?
        ##        f <- fitdistr(gw, "negative binomial", lower=c(0.001,0.01))
        f <- fitdistr(gw, "negative binomial")
        logp <- pnbinom(gw, size=f$estimate[['size']], mu=f$estimate[['mu']],
                        lower.tail=FALSE, log.p=TRUE)
        return(logp)
    } else {
        return(rep_len(0, length(gw)))
    }
}

geneNonPauseFitParam <- function(gw, avgcut=0.1, maxcut=4) {
    require(MASS)
    if (mean(gw)>avgcut & max(gw)>maxcut) {
        f <- fitdistr(gw, "negative binomial")
        f$estimate[['size']]
    } else
        1000
}

geneListNonPauseLogP <- function(gwl) {
    lapply(gwl, geneNonPauseLogP)
}



plotGeneFittedDistribution <- function(gw, ...) {
    require(MASS)
    f <- fitdistr(gw, "negative binomial")
#    t <- tabulate(gw)
#    plot(t, type='s')
    plot(density(gw),log='x', ...)
    curve(dnbinom(floor(x), size=f$estimate[['size']], mu=f$estimate[['mu']]),add=TRUE)
}

plotGene <- function(gw, ...) {
    plot(gw,type='h', ...)
}

pausingTable <- function(gwlp,gwlc, cds, trseq=polyatranscripts) {
    ss <- stack(gwlp)
    ss$gpos <- stack(lapply(gwlp, seq_along))$values
    names(ss) <- c('logP', 'gene', 'gpos')
    ss$count <- stack(gwlc)$values
    genelength <- sum(width(cds))
    genecounts <- sapply(gwlc,sum)
    ss$len <- genelength[ss$gene]
    ss$genecount <- genecounts[ss$gene]
    ss <- ss[order(ss$logP),]
    ss$logQ <- ss$logP-log(seq_len(nrow(ss)))+log(nrow(ss))
    ss$logP <- ss$logP*log10(exp(1))
    ss$logQ <- ss$logQ*log10(exp(1))
    ss$codon <- as.character(subseq(trseq[ss$gene], start=(ss$gpos-1)*3+1, width=3))
    ss$acodon <- as.character(subseq(trseq[ss$gene], start=(ss$gpos)*3+1, width=3))
    ss
}
