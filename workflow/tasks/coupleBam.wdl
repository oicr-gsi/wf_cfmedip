version 1.0

task coupleBam{
  input{
    String bamFiltered1
    String bamFiltered2
    String outputPath
    String fname
    Int threads
  }
  
  String bamBowtie2Subtracted=outputPath+"/"+fname+".bowtie2.subtracted.bam"
  String readsMagicblast=outputPath+"/read_ids_magic-blast.txt"
  String bamHybrid=outputPath+"/"+fname+".hybrid.bam"
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    
    samtools view -@ 8 ~{bamFiltered1} | awk '~{bracketOpen} print $1 ~{bracketClose}' > ~{readsMagicblast}
    
    java -jar /usr/lib/picard.jar FilterSamReads \
    I=~{bamFiltered2} \
    O=~{bamBowtie2Subtracted} \
    READ_LIST_FILe=~{readsMagicblast} \
    FILTER=excludeReadList
    
    samtools merge -@ ~{threads} ~{bamHybrid} ~{bamFiltered1} ~{bamBowtie2Subtracted}
    
  }
  
  output{
    File bamHybrid=bamHybrid
  }
}