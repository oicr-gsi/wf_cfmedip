library(docopt)

doc <- "Usage:
test.r --bamFile <FILE> --outputDir <DIR> --windowSize <SIZE>
runMedips.r --bamFile=FILE --outputDir=DIR --windowSize=SIZE

--bamFile FILE       Aligned, sorted, filtered reads (bam) [default: ]
--outputDir DIR      Path to output folder [default: ]
--windowSize SIZE    Size of genomic windows (bp) [default: ]
--help               show this help text"
opt <- docopt(doc)

if (!file.exists(opt$bamFile)){
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
ws=as.numeric(opt$windowSize)
paired=TRUE

# Disables the scientific notation to avoid powers in genomic coordinates (i.e. 1e+10)
options(scipen = 999)

chr.select=paste0("chr",c(1:22,"X","Y"))

message("MEDIPS: ws:",ws)
MeDIPset = MEDIPS.createSet(file = opt$bamFile,BSgenome = BSgenome, extend = extend, shift = shift, paired=paired, uniq = uniq, window_size = ws,chr.select=chr.select)

fname<-unlist(strsplit(basename(opt$bamFile),split="\\."))[1]
if(!dir.exists(opt$outputDir)){dir.create(opt$outputDir)}

pos.start<-c()
pos.end<-c()
pos.chr<-c()
for(i in 1:length(chr.select)){
  no.window<-ceiling(MeDIPset@chr_lengths[i]/ws)
  pos.chr<-c(pos.chr,rep(chr.select[i],no.window))
  pos.start<-c(pos.start,seq(from=1,to=((no.window-1)*ws)+1,by=ws))
  pos.end<-c(pos.end,seq(from=ws,to=no.window*ws,by=ws))
}
pos.chr<-data.frame(pos.chr)
pos.start<-data.frame(pos.start)
pos.end<-data.frame(pos.end)
genome.counts<-data.frame(MeDIPset@genome_count)

df.counts<-cbind(pos.chr,pos.start,pos.end,genome.counts)
file.counts<-paste0(opt$outputDir,"/MEDIPS_hg38_",fname,"_ws",ws,"_count.txt")
write.table(df.counts,file.counts,row.names=F,quote=F,col.names=F)
system(paste0("gzip -f ",file.counts))

file.rms<-paste0(opt$outputDir,"/MEDIPS_hg38_",fname,"_ws",ws,"_rms.txt")
#coupling set: maps CG densities across the genome
CS=MEDIPS.couplingVector(pattern="CG", refObj=MeDIPset)
df.rms<-NULL
tryCatch({
  #####MEDIPS.set(MeDIP=true,...) performs CpG density normalization (rms: relative methylation score)
  MEDIPSset.meth=MEDIPS.meth(MSet1 = MeDIPset, CSet = CS,MeDIP=TRUE)
  MEDIPS.rms<-MEDIPSset.meth[,paste0(gsub("-",".",basename(opt$bamFile)),".rms")]
  df.rms<-MEDIPSset.meth[,c('chr','start','stop',paste0(gsub("-",".",basename(opt$bamFile)),".rms"))]
},error=function(e){
  message("Error: MEDIPS CpG density normalization failed due to small number of reads")
})
if(is.null(df.rms)){df.rms<-("#Error: MEDIPS CpG density normalization failed due to small number of reads")}
write.table(df.rms,file.rms,row.names=F,quote=F,col.names=F)
system(paste0("gzip -f ",file.rms))

#CpG enrichment
#    Performance can be improved, it currently performs a full MEDIPS run. 
er <- MEDIPS.CpGenrich(file=opt$bamFile, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, chr.select=chr.select)
df.er<-data.frame(matrix(unlist(er), nrow=length(1)))
colnames(df.er)<-names(er)
df.er<-cbind(sample=fname,df.er)
file.er<-paste0(opt$outputDir,"/MEDIPS_hg38_",fname,"_CpGenrich.txt")
write.table(df.er,file.er,row.names=F,quote=F,col.names=F)

