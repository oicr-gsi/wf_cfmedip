version 1.0

task extractUMI {
  input{
    Boolean useUMI
    File R1
    File R2
    String outputPath
    String fname
    String patternUMI
    String patternUMI2
    Int newReadLen
  }
  
  command{
    mkdir -p ~{outputPath}
    mkdir -p ~{outputPath}/metrics
    
    if [[ ~{useUMI} == true ]];then  
    umi_tools extract --extract-method=string \
    --bc-pattern=~{patternUMI} --bc-pattern2=~{patternUMI2} \
    -I ~{R1} \
    --read2-in=~{R2} \
    -S ~{outputPath}/~{fname}.R1.fq.gz \
    --read2-out=~{outputPath}/~{fname}.R2.fq.gz \
    -L ~{outputPath}/metrics/extractUMI.log
    
    else
      cp ~{R1} ~{outputPath}/~{fname}.R1.fq.gz
    cp ~{R2} ~{outputPath}/~{fname}.R2.fq.gz
    
    fi
    
    if [[ ~{newReadLen} != -1 ]];then
    fastp \
    --max_len1 ~{newReadLen} \
    --max_len2 ~{newReadLen} \
    -i ~{outputPath}/~{fname}.R1.fq.gz \
    -I ~{outputPath}/~{fname}.R2.fq.gz \
    -o ~{outputPath}/~{fname}.R1_trim.fq.gz \
    -O ~{outputPath}/~{fname}.R2_trim.fq.gz \
    --disable_adapter_trimming \
    --disable_trim_poly_g \
    --disable_quality_filtering \
    --json ~{outputPath}/metrics/fastp_trim.json \
    --html /dev/null
    
    mv -f ~{outputPath}/~{fname}.R1_trim.fq.gz ~{outputPath}/~{fname}.R1.fq.gz
    mv -f ~{outputPath}/~{fname}.R2_trim.fq.gz ~{outputPath}/~{fname}.R2.fq.gz
    fi
  }
  
  output {
    File extrR1=outputPath+"/"+fname+".R1.fq.gz"
    File extrR2=outputPath+"/"+fname+".R2.fq.gz"
  }
}