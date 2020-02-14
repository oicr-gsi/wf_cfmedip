version 1.0

task removeDuplicates{
  input{
    Boolean useUMI
    String bamFiltered
    String fname
    String outputPath
    Int threads
  }
  
  String bamDedup=outputPath+"/"+fname+".dedup.bam"
  
  #index bamFilter required by umi_tools; index filtered.dedup.bam for alignment visualization
  
  command{   
    if [[ ~{useUMI} == true ]];then
    samtools index -@ ~{threads} ~{bamFiltered}
    
    umi_tools dedup --paired \
    -I ~{bamFiltered} \
    -S ~{bamDedup} \
    --output-stats=~{outputPath}/metrics/removeDuplicates-UMI_tools 
    
    else
      java -jar /usr/lib/picard.jar MarkDuplicates \
    I=~{bamFiltered} \
    O=~{bamDedup} \
    M=~{outputPath}/metrics/removeDuplicates-statsPicard.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=true
    
    fi
    
    samtools index -@ ~{threads} ~{bamDedup}
  }
  
  output{
    File bamDedup=bamDedup
  }
}