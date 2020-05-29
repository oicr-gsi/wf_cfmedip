version 1.0

task alignBwa{
  input{
    String extrR1
    String extrR2
    String indexBwa
    String fname
    String outputPath
    Int threads
    
    #holds bwa mem until magicblast finishes and threads are available
    File? processHold
  }
  
  String aligner="bwa"
  
  command{
    bwa mem -t ~{threads} \
    ~{indexBwa} \
    ~{extrR1} \
    ~{extrR2} \
    > ~{outputPath}/~{fname}.bwa.sam
    samtools view -@ ~{threads} -b ~{outputPath}/~{fname}.~{aligner}.sam | samtools sort -o ~{outputPath}/~{fname}.~{aligner}.bam
    rm ~{outputPath}/~{fname}.~{aligner}.sam

  }
  
  output{
    File bamAligned=outputPath+"/"+fname+"."+aligner+".bam"
  }
}
