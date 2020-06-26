version 1.0

task runMedestrand{
  input{
    File bamDedup
    String fname
    String outputPath
    Int windowSize
    Boolean useMedestrand
  }
  
  String outMedips=outputPath+'/runMedips'
  
  command{
  if [[ ~{useMedestrand} == true ]];then
    mkdir -p ~{outMedips}
    r /workflow/runMedips/runMedestrand.r \
    --bamFile ~{bamDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
    fi
  }
  output{
    Array[File] optional_medestrandWig = glob(outMedips+'/MeDESTrand_hg38_'+fname+'_ws'+windowSize+'_wig.bed.gz')
  }
}
