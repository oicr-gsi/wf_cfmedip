version 1.0

task runMedips{
  input{
    File bamDedup
    File ROIFile
    String fname
    String outputPath
    Int windowSize
  }
  
  String outMedips=outputPath+'/runMedips'
  String ROIName=sub(basename(ROIFile),"(\.bed)?", "")
  
  command{
    mkdir -p ~{outMedips}
    
    r /workflow/runMedips/runMedips.r \
    --bamFile ~{bamDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
    
    r /workflow/runMedips/runMedipsROI.r \
    --bamFile ~{bamDedup} \
    --outputDir ~{outMedips} \
    --ROIFile ~{ROIFIle}
    
  }
  output{
    File medipsCount=outMedips+'/MEDIPS_hg38_'+fname+'_ws'+windowSize+'_count.txt.gz'
    File medipsRms=outMedips+'/MEDIPS_hg38_'+fname+'_ws'+windowSize+'_rms.txt.gz'
    File medipsROI=outMedips+'/'+fname+'_'+ROIName+'_CPM.RData'
  }
}




