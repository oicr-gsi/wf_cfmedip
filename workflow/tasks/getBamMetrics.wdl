version 1.0

task getBamMetrics{
  input{
    String fname
    String outputPath
    String bamDedup
    String fasta
  }
  
  String picardOut=outputPath+'/picard'
  
  command{
    mkdir -p ~{picardOut}
    
    java -jar /usr/lib/picard.jar CollectMultipleMetrics \
    PROGRAM=null \
    PROGRAM=CollectGcBiasMetrics \
    PROGRAM=CollectInsertSizeMetrics \
    R=~{fasta} \
    I=~{bamDedup} \
    O=~{picardOut}/picard \
    VALIDATION_STRINGENCY=SILENT    
  }
  
  output{
    File picardInsertSizeMetrics=picardOut+'/picard.insert_size_metrics'
    File picardGcBiasMetrics=picardOut+'/picard.gc_bias.summary_metrics'
  }
  
}