#Adapted from: https://github.com/AdaZED/medestrand/blob/master/medestrand_counts_generator.R
library(optparse)
msg.usage="Rscript %prog --inputDir input/dir --outputDir output/dir --bamFile sample.bam"
msg.description="Extracts genomic window counts in cfMeDIP"
option_list = list(
  make_option("--bamFile", type="character", default=NULL, help="REQUIRED\tAligned, sorted, filtered reads (bam)", metavar="character"),
  make_option("--outputDir", type="character", default=NULL, help="REQUIRED\tPath to output folder", metavar="character"),
  make_option("--windowSize", type="integer", default=NULL, help="REQUIRED\tSize of genomic windows (bp)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,description=msg.description,usage = msg.usage);
opt = parse_args(opt_parser);

if (is.null(opt$bamFile) | is.null(opt$outputDir | is.null(opt$windowSize))){
  print_help(opt_parser)
  stop("ERROR: parameters --bamFile, --outputDir and --windowSize must be provided", call.=FALSE)
}
if (!file.exists(paste0(opt$inputDir,"/",opt$bamFile))){
  print_help(opt_parser)
  stop(paste0("ERROR: bam file not found ",opt$bamFile), call.=FALSE)
}
if (!file.exists(opt$outputDir)){
  dir.create(opt$outputDir)
}

library(MeDEStrand)
library("BSgenome.Hsapiens.UCSC.hg38")
library(GenomicRanges)
library(gsubfn)
args=(commandArgs(TRUE))

# Retrieve the user parameters
sample <- opt$bamFile
output <- opt$outputDir
ws <- opt$windowSize
paired <- TRUE

#  Adapted from: https://github.com/jxu1234/MeDEStrand/blob/master/R/MeDEStrand.createSet.R
#  The original function  uses hardcoded hg19; here, we switch to hg38
MeDEStrand.binMethyl_hg38 <- function(MSetInput=NULL, CSet=NULL, ccObj=NULL, Granges = FALSE){
  for ( i in 1:2 ) {
    if(is.list(MSetInput)){
      MSet=MSetInput[[i]]
    }
    signal =  genome_count(MSet)
    coupling = genome_CF(CSet)
    ccObj = MeDEStrand.calibrationCurve(MSet=MSet, CSet=CSet, input=F)
    index.max = which(ccObj$mean_signal== max( ccObj$mean_signal[1:ccObj$max_index] ) )
    MS = ccObj$mean_signal[1:index.max]
    CF = ccObj$coupling_level[1:index.max]
    model.data = data.frame( model.MS =  MS/max( MS), model.CF = CF  )
    logistic.fit = glm(  model.MS ~ model.CF , family=binomial(logit) , data = model.data)
    if ( i == 1) { cat("Estimating and correcting CG bias for reads mapped to the DNA positive strand...\n") }
    if ( i == 2) { cat("Estimating and correcting CG bias for reads mapped to the DNA negative strand...\n") }
    estim=numeric(length(ccObj$mean_signal))  # all 0's
    low_range=1:index.max
    estim[low_range]=ccObj$mean_signal[low_range]
    high_range = ( length(low_range)+1 ):length(estim)
    y.predict = predict( logistic.fit, data.frame( model.CF = ccObj$coupling_level[high_range] ), type ="response"  )*ccObj$mean_signal[ccObj$max_index]
    estim[high_range] = y.predict
    signal=signal/estim[coupling+1]
    signal[coupling==0]=0
    signal = log2(signal)
    signal[is.na(signal)] = 0
    minsignal=min(signal[signal!=-Inf])
    signal[signal!=-Inf]=signal[signal!=-Inf]+abs(minsignal)
    maxsignal = quantile(signal[signal!=Inf], 0.9995  )
    signal[signal!=Inf & signal>maxsignal]=maxsignal
    signal=round((signal/maxsignal), digits=2)
    signal[signal==-Inf | signal ==Inf]=0
    if ( i == 1) {  pos.signal = signal }
    if ( i == 2) { neg.signal = signal }
  }
  merged.signal = (pos.signal+neg.signal)/2
  if( !Granges ) {
    return(  merged.signal)}else{
      chr.select = MSet@chr_names
      window_size = window_size(MSet)
      chr_lengths=unname( seqlengths(BSgenome.Hsapiens.UCSC.hg38)[ seqnames(BSgenome.Hsapiens.UCSC.hg38@seqinfo)%in%chr.select ] )
      no_chr_windows = ceiling(chr_lengths/window_size)
      supersize_chr = cumsum(no_chr_windows)
      chromosomes=chr.select
      all.Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)
      all.Granges.genomeVec$CF = CS@genome_CF
      all.Granges.genomeVec$binMethyl= merged.signal
      return( all.Granges.genomeVec )
    }
}


# Disable the scientific notation in R (to avoid powers (1e+10 for instance) to be written in output files)
options(scipen = 999)

# Set global variables for importing short reads. For details, in R console, type "?MeDEStrand.createSet"
BSgenome="BSgenome.Hsapiens.UCSC.hg38"
uniq = 1
extend = 200
shift = 0
chr.select = paste0('chr', c(1:22) )

# Create a MeDIP set
MeDIP_seq = MeDEStrand.createSet(file=paste0(opt$inputDir,"/",opt$bamFile), BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr.select, paired=paired)

#  Count CpG pattern in the bins
CS = MeDEStrand.countCG(pattern="CG", refObj=MeDIP_seq)

# Infer genome-wide absolute methylation levels:
#result.methylation = MeDEStrand.binMethyl(MSetInput = MeDIP_seq, CSet = CS, Granges = TRUE)
result.methylation = MeDEStrand.binMethyl_hg38(MSetInput = MeDIP_seq, CSet = CS, Granges = TRUE)

# Create a dataframe from the previous GRanges object.
# Warning: GRanges and UCSC BED files use different conventions for the genomic coordinates
# GRanges use 1-based intervals (chr1:2-8 means the 2nd till and including the 8th base of chr1, i.e. a range of length of 7 bases)
# UCSC bed-files use 0-based coordinates (chr1:2-8 means the 3rd base till and including the 8th base, i.e. a range of length of 6 bases)

# Dataframe for generating a bed file used to generate then a big matrix from all the samples
df_for_matrix <- data.frame(seqnames=seqnames(result.methylation),
                            starts=start(result.methylation)-1,
                            ends=end(result.methylation),
                            scores=elementMetadata(result.methylation)$binMethyl)

# Dataframe for generating a bed file used to generate then a wig file
df_for_wig <- data.frame(seqnames=seqnames(result.methylation),
                         starts=start(result.methylation)-1,
                         ends=end(result.methylation),
                         names=c(rep(".", length(result.methylation))),
                         scores=elementMetadata(result.methylation)$binMethyl,
                         strands=strand(result.methylation))

# Export in a bed file for the final matrix
extformatrix <- "_matrix.bed"
bed_matrix_output <- fn$identity("$output$sample$extformatrix")
write.table(df_for_matrix, file=bed_matrix_output, quote=F, sep="\t", row.names=F, col.names=F)

# Export in a bed file for the wig files
extforwig <- "_wig.bed"
bed_wig_output <- fn$identity("$output$sample$extforwig")
write.table(df_for_wig, file=bed_wig_output, quote=F, sep="\t", row.names=F, col.names=F)

# Export in a wig file
wig <- ".wig"
bed_loaded <- import(con=bed_wig_output, format="bed")
wig_output <- fn$identity("$output$sample$wig")
export.wig(object=bed_loaded, con=wig_output)
