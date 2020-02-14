version 1.0

task alignBowtie2{
  input{
    String extrR1
    String extrR2
    String indexBowtie2
    String fname
    String outputPath
    Int threads
    
    #holds bowtie2 until magicblast finishes and threads are available
    File? processHold
  }
  
  String aligner="bowtie2"
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    
    bowtie2 -p ~{threads} -x ~{indexBowtie2} \
    -1 ~{extrR1} \
    -2 ~{extrR2} \
    -S ~{outputPath}/~{fname}.bowtie2.sam
    samtools view -@ ~{threads} -b ~{outputPath}/~{fname}.~{aligner}.sam | samtools sort -o ~{outputPath}/~{fname}.~{aligner}.bam
    rm ~{outputPath}/~{fname}.~{aligner}.sam
    
  }
  
  output{
    File bamAligned=outputPath+"/"+fname+"."+aligner+".bam"
  }
}