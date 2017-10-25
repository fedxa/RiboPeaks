# Example peak search analysis
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Matrix)
library(MASS)
source('peaksearch.R')

## Name of the alignment BAM file
bams <- "alignment.bam"

## Read the genome
##   Genome can be obtained from https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
##   We use the fasta given in in the saccharomyces_cerevisiae_R64-1-1_20110208.gff
yg <- FaFile("reference/yeast.fa")
## Read the gene annotations
txdb <- makeTxDbFromGFF("reference/yeast_verified_genes.gff")

cds <- cdsBy(txdb, by="gene")

## Extract the transcripts for genes
## Do not use!!!! does not rearange multiple exons on '-' strand
### transcripts.bad <- extractTranscriptSeqs(yg,cds)
## Correct, but a bit slow function is
transcripts <- properExtractTranscripts(yg, cds)

## Augment the transcripts by AAA on each end
##  allows to safely look at codons before or after teh pausing codon
polyatranscripts <- polyaTranscripts(transcripts)

## Get the list of P-sites from RPF fragments of slected lengths
## with given offset into the fragment
##  Check https://github.com/fedxa/RPF_pipeline
##  for semiautomatic way
psites <- readPsites(bams, c(27,28,29,30), c(11,12,12,12))

wl <- sitesToWigVector(psites)
geneWigs <- extractGeneListWigs(cds, wl)
geneWigsC <- codonWigs(geneWigs)
geneWigsLogP <- geneListNonPauseLogP(geneWigsC)

## Produces the table of candidate pausing sites
pt <- pausingTable(geneWigsLogP, geneWigsC, cds)


## Let us convert this to the genomic coordinates
ptWigByGene <- with(pt,GRanges(seqnames=gene, IRanges(start=(gpos-1)*3+1,width=3)))
ptWig <- mapFromTranscripts(ptWigByGene,cds)
mcols(ptWig) <- pt[ptWig$xHits,c("gene","logP","count","len","genecount","logQ","codon","gpos","acodon")]
ptWig$codon<-factor(as.character(ptWig$codon))
ptWig$acodon<-factor(as.character(ptWig$acodon))
ptWig$ratio<-ptWig$count*ptWig$len/ptWig$genecount

## Select only reasonable expressed genes -- RPKM>0.1
ptWigExp <- ptWig[ptWig$genecount/ptWig$len*3 > 0.1]
ptWigExp <- ptWigExp[order(ptWigExp$ratio,decreasing=TRUE)]


# Save the final table as a CSV file
write.csv(ptWigExp, "pausing.csv")
