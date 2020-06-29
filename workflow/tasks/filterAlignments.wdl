version 1.0

task filterAlignments{
  input{
    File bamAligned
    String fname
    String outputPath
    String aligner
    Int maxMismatch
    Int threads
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{

    samtools view -@ ~{bamAligned} \
    | awk 'sqrt($9*$9)>119 && sqrt($9*$9)<501' \
    | awk '~{bracketOpen}print $1~{bracketClose}' \
    > ~{outputPath}/~{fname}.~{aligner}.filter1.mapped_proper_pair.txt
    
    java -jar /usr/lib/picard.jar FilterSamReads \
    I=~{bamAligned} \
    O=~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    READ_LIST_FILE=~{outputPath}/~{fname}.~{aligner}.filter1.mapped_proper_pair.txt \
    FILTER=includeReadList
    
    samtools view -@ ~{threads} ~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    | awk '~{bracketOpen}read=$0;sub(/.*NM:i:/,X,$0);sub(/\t.*/,X,$0);if(int($0)>~{maxMismatch})~{bracketOpen}print read~{bracketClose}~{bracketClose}' \
    | awk '~{bracketOpen}print $1~{bracketClose}' \
    > ~{outputPath}/~{fname}.~{aligner}.filter2.high_mismatch.txt
    
    java -jar /usr/lib/picard.jar FilterSamReads \
    I=~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    O=~{outputPath}/~{fname}.~{aligner}.filter3.bam \
    READ_LIST_FILe=~{outputPath}/~{fname}.~{aligner}.filter2.high_mismatch.txt \
    FILTER=excludeReadList
    
  }
  
  output{
    File bamFiltered=outputPath+"/"+fname+"."+aligner+".filter3.bam"
  }
}
