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

# Disables the scientific notation to avoid powers in genomic coordinates (i.e. 1e+10)
options(scipen = 999)

chr.select=paste0("chr",c(1:22,"X","Y"))

message("MEDIPS: ws:",ws,".txt")
MeDIPset = MEDIPS.createSet(file = opt$bamFile,BSgenome = BSgenome, extend = extend, shift = shift, paired=paired, uniq = uniq, window_size = ws,chr.select=chr.select)

#coupling set: maps CG densities across the genome
CS=MEDIPS.couplingVector(pattern="CG", refObj=MeDIPset)

#####MEDIPS.set(MeDIP=true,...) performs CpG density normalization (rms: relative methylation score)
MEDIPSset.meth=MEDIPS.meth(MSet1 = MeDIPset, CSet = CS,MeDIP=TRUE)

fname<-unlist(strsplit(basename(opt$bamFile),split="\\."))[1]
if(!dir.exists(paste0(opt$outputDir,"/runMedips"))){dir.create(paste0(opt$outputDir,"/runMedips"))}

df.counts<-cbind(MEDIPSset.meth[,c('chr','start','stop')],MeDIPset@genome_count)
file.counts<-paste0(opt$outputDir,"/runMedips/MEDIPS_hg38_",fname,"_ws",ws,"_count.txt")
write.table(df.counts,file.counts,row.names=F,quote=F,col.names=F)
system(paste0("gzip ",file.counts))

MEDIPS.rms<-MEDIPSset.meth[,paste0(gsub("-",".",basename(opt$bamFile)),".rms")]
df.rms<-MEDIPSset.meth[,c('chr','start','stop',paste0(gsub("-",".",basename(opt$bamFile)),".rms"))]
file.rms<-paste0(opt$outputDir,"/runMedips/MEDIPS_hg38_",fname,"_ws",ws,"_rms.txt")
write.table(df.rms,file.rms,row.names=F,quote=F,col.names=F)
system(paste0("gzip ",file.rms))










