library(optparse)
library(parallel)
library(reshape2)


msg.usage="Rscript %prog --study study_example [additional options]"
msg.description="Wrapper for the minfi package to extract and normalize data from idat files. A SampleSheet file must be present"
option_list = list(
  make_option("-r1", type="character", default=NULL, help="REQUIRED\tFile read1.fastq.gz", metavar="character"),
  make_option("-r2", type="character", default=NULL, help="REQUIRED\tFile read2.fastq.gz", metavar="character"),
  make_option("-index", type="character", default=NULL, help="REQUIRED\tBowtie2 index", metavar="character"),
  make_option("-crop", type="character", default=NULL,help="default=NULL\tOPTIONAL\tRead cropping? Takes Trimmomatic parameters as input", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,description=msg.description,usage = msg.usage);
opt = parse_args(opt_parser);

if (is.null(opt$r1) | is.null(opt$r2)){
  print_help(opt_parser)
  stop("ERROR: parameter -r1 and -r2 must be provided", call.=FALSE)
}
if (!file.exists(opt$r1)){
  print_help(opt_parser)
  stop(paste0("ERROR: file not found -r1 ",opt$r1), call.=FALSE)
}
if (!file.exists(opt$r2)){
  print_help(opt_parser)
  stop(paste0("ERROR: file not found -r2 ",opt$r1), call.=FALSE)
}
if (is.null(opt$index)){
  print_help(opt_parser)
  stop("ERROR: parameter -index must be provided", call.=FALSE)
}

fname<-system(paste0("basename ",opt," | sed 's/.R1.*//'"),intern=TRUE)
bt2.index<-opt$index

#step 0: create directories
dirs<-c("1_trim","2_align_qc","3_preprocess","4_picard","5_consensuscruncher","6_medips")
sapply(dirs,function(x) system(paste0("mkdir -p ",x)))


#step 1: trim
if(!is.null(opt$crop)){
  cmd<-paste0("/usr/bin/TrimmomaticPE ",opt$r1," ",opt$r2," 1_trim/",fname,".R1.fq.gz 1_trim/",fname,".R1_unpaired.fq.gz 1_trim/",fname,".R2.fq.gz 1_trim/",fname,".R2_unpaired.fq.gz ",opt$crop)
  system(cmd)
}else{
  system(paste0("ln -s ",opt$r1," 1_trim/",fname,"R1.fq.gz"))
  system(paste0("ln -s ",opt$r2," 1_trim/",fname,"R2.fq.gz"))
}

#step 2: alignment
cmd<-paste0("bowtie2 -p 8 -x ",bt2.index," -1 1_trim/",fname,".R1.fq.gz -2 1_trim/",fname,".R2.fq.gz -S 2_align_qc/",fname,".sam")
system(cmd)

#step 3: preprocessing
#samtools view -b (output bam) -F12 (remove unpaired reads)
cmd<-paste0("samtools view -b -F12 2_align_qc/",fname,".sam | samtools sort -o 3_preprocess/",fname,".sorted.bam")
system(cmd)


cmd<-paste0("java -jar /usr/bin/picard.jar MarkDuplicates -I 3_preprocess/",fname,".sorted.bam -O 3_preprocess/",
            fname,".sorted.dedup.bam"," -M 3_preprocess/",fname,".sorted.dedup.metrics -ASSUME_SORTED true -VALIDATION_STRINGENCY SILENT -REMOVE_DUPLICATES true")
system(cmd)



#step 4: get some alignment metrics
cmd<-paste0("java -jar /usr/bin/picard.jar CollectMultipleMetrics R=data/Bowtie2/",bt2.index,".fa I=3_preprocess/",fname,".sorted.dedup.bam",
            " O=4_picard/",fname," VALIDATION_STRINGENCY=SILENT")
system(cmd)


cmd<-paste0("java -jar /usr/bin/picard.jar CollectGcBiasMetrics R=data/Bowtie2/",bt2.index,".fa I=3_preprocess/",fname,".sorted.dedup.bam",
            " O=4_picard/",fname,".gc_bias_metrics.txt S=4_picard/",fname,".summary_gc_bias_metrics.txt CHART=4_picard/",fname,".gc_bias_metrics.pdf")
system(cmd)



#step 4: MEDIPS counts
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=0 ####WARNING: default settings normalize the data, must be set to 0 to disable this transformation
extend=0 #relevant for single-end: https://support.bioconductor.org/p/81098/
shift=0
list.ws=c(100,200,300)
paired=TRUE
chr.select=paste0("chr",c(1:22,"X","Y"))
for(ws in list.ws){
  message("MEDIPS: ",fname,".ws",ws,".txt")
  MeDIP.set = MEDIPS.createSet(file = paste0("3_preprocess","/",fname,".sorted.dedup.bam"),BSgenome = BSgenome, extend = extend, shift = shift, paired=paired, uniq = uniq, window_size = ws,chr.select=chr.select)
  write.table(data.frame(MeDIP.set@genome_count),paste0("6_medips/",fname,".ws",ws,".txt"),row.names=F,quote=F,col.names=F)
}


#step 6: normalize by CpG
load(file=paste0("../data/window_name_island_hs37_ws",ws,".RData"))
is.island<-rep(FALSE,length(name.island))
is.island[name.island!="none"]<-TRUE
sum(is.island)
load(file=("../data/CpG_island_EPIC.RData"))

fname<-"CMP-01-03-cfDNA-03.txt"
message(fname,"...")

counts<-read.table(paste0("6_medips/",fname,".ws",ws,".txt"),header=F)[,1]
counts.sub<-counts[is.island]
name.island.sub<-name.island[is.island]


cl <- makeCluster(4)
clusterExport(cl, list("counts.sub","name.island.sub"))
island.count<-parSapply(cl, island, function(x) {
  counts.curr<-counts.sub[name.island.sub==x]
  return(sum(counts.curr))
})
stopCluster(cl)


chunks<-colsplit(colsplit(island,":",letters[1:2])$b,"-",letters[2:3])
island.start<-chunks$b
island.end<-chunks$c

#https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
TPM.RPK<-island.count/((island.end-island.start)*0.001)
TPM<-TPM.RPK/(sum(TPM.RPK)/1e6)
#we consider these FPKM, as MEDIPS only counts one hit per read pair, and single-mapped reads were removed
FPKM<-(island.count/(sum(counts)/1e6))/((island.end-island.start)*0.001)


df<-data.frame(island=island,count=island.count,FPKM=FPKM,TPM=TPM)
write.csv(df,paste0("6_medips/",fname,".ws",ws,".CpG_norm.txt"),row.names = FALSE,quote=F)







