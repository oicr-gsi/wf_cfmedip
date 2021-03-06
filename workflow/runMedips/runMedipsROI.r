library(docopt)

doc <- "Usage:
runMedipsROI.r --bamFile <FILE> --ROIFile <FILE> --outputDir <DIR>
runMedipsROI.r --bamFile=FILE --ROIFile=FILE --outputDir=DIR
--bamFile FILE       Aligned, sorted, filtered reads (bam) [default: ]
--ROIFile FILE       Regions of interest, with columns: chr,start,end,name [default: ]
--outputDir DIR      Path to output folder [default: ]
--help               Show this help text"
opt <- docopt(doc)

if (!file.exists(opt$bamFile)){
  stop(paste0("ERROR: bam file not found ",opt$bamFile), call.=FALSE)
}
if (!file.exists(opt$ROIFile)){
  stop(paste0("ERROR: ROI file not found ",opt$ROIFile), call.=FALSE)
}
if (!file.exists(opt$outputDir)){
  dir.create(opt$outputDir)
}
sample.name<-unlist(strsplit(basename(opt$bamFile),split="\\."))[1]
ROI.name<-unlist(strsplit(basename(opt$ROIFile),split="\\."))[1]

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)
library(doParallel)

BSgenome="BSgenome.Hsapiens.UCSC.hg38"
uniq=0 ####WARNING: default settings normalize the data, must be set to 0 to disable this transformation
extend=0 #relevant for single-end: https://support.bioconductor.org/p/81098/
shift=0
paired=TRUE

# Disables the scientific notation to avoid powers in genomic coordinates (i.e. 1e+10)
options(scipen = 999)

chop<-function(df,block.size=10000){
  if(!is.data.frame(df)){df<-as.data.frame(df)}
  if(block.size>=nrow(df)){
    block.size<-nrow(df) %/% 2 #integer division
  }
  no.blocks<-nrow(df) %/% block.size #integer division
  block.map<-NULL
  c<-1
  for(i in 1:no.blocks){
    block.map<-c(block.map,rep(c,block.size))
    c<-c+1
  }
  block.map<-c(block.map,rep(c,nrow(df)-(block.size*no.blocks)))
  chunks<-split(df, block.map)
  return(chunks)
}

chr.select=paste0("chr",c(1:22,"X","Y"))

ROI<-read.table(opt$ROIFile,header=TRUE,stringsAsFactors = FALSE)

#does not work inside docker, due to ssh calls ...
#cl <- makeCluster(opt$threads)
#clusterEvalQ(cl, library("MEDIPS"))
#doParallel::registerDoParallel(cl)
#ROI.counts<-foreach(ROI.block=chop(ROI,block.size = 10000),.combine='c') %dopar% {
#  library(BSgenome.Hsapiens.UCSC.hg38) # problem to register this library in the cluster, needs to be loaded directly by the nodes
#  MEDIPS.createROIset(file=opt$bamFile, ROI=ROI.block, extend=extend,
#                      shift=shift, bn=1, BSgenome="BSgenome.Hsapiens.UCSC.hg38", uniq=uniq, 
#                      chr.select=chr.select[chr.select %in% unique(ROI.block$chr)], 
#                      paired = paired, sample_name=NULL, 
#                      isSecondaryAlignment=FALSE, simpleCigar=TRUE)@genome_count
#}
#stopCluster(cl)



ROI.counts<-foreach(ROI.block=chop(ROI,block.size = 10000),.combine='c') %do% {
  MEDIPS.createROIset(file=opt$bamFile, ROI=ROI.block, extend=extend,
                      shift=shift, bn=1, BSgenome="BSgenome.Hsapiens.UCSC.hg38", uniq=uniq, 
                      chr.select=chr.select[chr.select %in% unique(ROI.block$chr)], 
                      paired = paired, sample_name=NULL, 
                      isSecondaryAlignment=FALSE, simpleCigar=TRUE)@genome_count
}





RPK<-ROI.counts/((ROI$end-ROI$start)*0.001)
CPM<-RPK/(sum(RPK)*0.000001)

ROI.signal<-data.frame(name=ROI$name,count=ROI.counts,CPM=CPM)
row.names(ROI.signal)<-NULL
fname<-paste0(opt$outputDir,"/",sample.name,"_",ROI.name,"_CPM.RData")    
save(ROI.signal,file=fname,compress=TRUE)



