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

if (is.null(opt$bamFile) | is.null(opt$outputDir) | is.null(opt$windowSize)){
  print_help(opt_parser)
  stop("ERROR: parameters --bamFile, --outputDir and --windowSize must be provided", call.=FALSE)
}
if (!file.exists(opt$bamFile)){
  print_help(opt_parser)
  stop(paste0("ERROR: bam file not found ",opt$bamFile), call.=FALSE)
}
if (!file.exists(opt$outputDir)){
  dir.create(opt$outputDir)
}

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

BSgenome="BSgenome.Hsapiens.UCSC.hg38"
uniq=0 ####WARNING: default settings normalize the data, must be set to 0 to disable this transformation
extend=0 #relevant for single-end: https://support.bioconductor.org/p/81098/
shift=0
ws=opt$windowSize
paired=TRUE

chr.select=paste0("chr",c(1:22,"X","Y"))

message("MEDIPS: ws:",ws,".txt")
MeDIP.set = MEDIPS.createSet(file = opt$bamFile,BSgenome = BSgenome, extend = extend, shift = shift, paired=paired, uniq = uniq, window_size = ws,chr.select=chr.select)
write.table(data.frame(MeDIP.set@genome_count),paste0(opt$outputDir,"/MEDIPS_hg38_ws",ws,"_count.txt"),row.names=F,quote=F,col.names=F)

#coupling set: maps CG densities across the genome
CS=MEDIPS.couplingVector(pattern="CG", refObj=MeDIP.set)

#####MEDIPS.set(MeDIP=true,...) performs CpG density normalization (rms: relative methylation score)
MEDIPS.set.rms=MEDIPS.meth(MSet1 = MeDIP.set, CSet = CS,MeDIP=TRUE)[paste0(gsub("-",".",basename(opt$bamFile)),".rms")]

write.table(MEDIPS.set.rms,paste0(opt$outputDir,"/MEDIPS_hg38_ws",ws,"_rms.txt"),row.names=F,quote=F,col.names=F)




