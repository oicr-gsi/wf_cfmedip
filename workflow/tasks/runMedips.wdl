version 1.0

task runMedips{
  input{
    File bamDedup
    String fname
    String outputPath
    Int windowSize
  }
  
  String outMedips=outputPath+'/runMedips'
  
  command{
    mkdir -p ~{outMedips}
    
    r /workflow/runMedips/runMedips.r \
    --bamFile ~{bamDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
    
    r /workflow/runMedips/runMedestrand.r \
    --bamFile ~{bamDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
  }
  output{
    File medipsCount=outMedips+'/MEDIPS_hg38_'+fname+'_ws'+windowSize+'_count.txt.gz'
    File medipsRms=outMedips+'/MEDIPS_hg38_'+fname+'_ws'+windowSize+'_rms.txt.gz'
    File medestrandWig=outMedips+'/MeDESTrand_hg38_'+fname+'_ws'+windowSize+'_wig.bed.gz'
  }
}
