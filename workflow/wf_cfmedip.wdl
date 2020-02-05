version 1.0

workflow wf_cfmedip {
  input{
    String outputPath
    String aligner
    String index
    File R1
    File R2
    File fasta
    
    String? sampleName
    String patternUMI = "NNNNN"
    String patternUMI2 = "NNNNN"
    String seqMeth = "F19K16"
    String seqUmeth = "F24B22"
    Boolean useUMI = true
    Int windowSize = 200
    Int threads = 4
    Int newReadLen = -1
  }
  
  String fname=if defined(sampleName) then select_first([sampleName,""]) else sub(basename(R1),"(\.fq)?(\.fastq)?(\.gz)?", "")
  
  call extractUMI {
    input:
      useUMI=useUMI,
      R1=R1,
      R2=R2,
      aligner=aligner,
      index=index,
      outputPath=outputPath,
      fname=fname,
      patternUMI=patternUMI,
      patternUMI2=patternUMI2,
      newReadLen=newReadLen
  }
  
  call alignReads {
    input:
      extrR1=extractUMI.extrR1,
      extrR2=extractUMI.extrR2,
      outputPath=outputPath,
      fname=fname,
      index=index,
      aligner=aligner,
      threads=threads
  }
  
  call filterBadAlignments{
    input:
      bamFile=alignReads.bamFile,
      outputPath=outputPath,
      fname=fname,
      aligner=aligner,
      threads=threads
  }
  
  call removeDuplicates{
    input:
      useUMI=useUMI,
      bamFilter=filterBadAlignments.bamFilter,
      outputPath=outputPath,
      fname=fname,
      aligner=aligner,
      threads=threads
  }
  
  call getBamMetrics{
    input:
      fname=fname,
      outputPath=outputPath,
      aligner=aligner,
      bamFilterDedup=removeDuplicates.bamFilterDedup,
      fasta=fasta
  }
  
  call parseMethControl{
    input:
      bamFilterDedup=removeDuplicates.bamFilterDedup,
      outputPath=outputPath,
      fname=fname,
      aligner=aligner,
      seqMeth=seqMeth,
      seqUmeth=seqUmeth,
      threads=threads
  }
  
  call runMedips{
    input:
      bamFilterDedup=removeDuplicates.bamFilterDedup,
      fname=fname,
      outputPath=outputPath,
      aligner=aligner,
      windowSize=windowSize
  }
  
}

task extractUMI {
  input{
    Boolean useUMI
    File R1
    File R2
    String aligner
    String index
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

task alignReads{
  input{
    String extrR1
    String extrR2
    String index
    String fname
    String outputPath
    String aligner
    Int threads
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    if [ "~{aligner}" == "bowtie2" ];then
    bowtie2 -p ~{threads} -x ~{index} \
    -1 ~{extrR1} \
    -2 ~{extrR2} \
    -S ~{outputPath}/~{fname}.bowtie2.sam
    samtools view -@ ~{threads} -b ~{outputPath}/~{fname}.~{aligner}.sam | samtools sort -o ~{outputPath}/~{fname}.~{aligner}.bam
    rm ~{outputPath}/~{fname}.~{aligner}.sam
    fi
       
    if [ "~{aligner}" == "bwa" ];then
    bwa mem -t ~{threads} \
    ~{index} \
    ~{extrR1} \
    ~{extrR2} \
    > ~{outputPath}/~{fname}.bwa.sam
    samtools view -@ ~{threads} -b ~{outputPath}/~{fname}.~{aligner}.sam | samtools sort -o ~{outputPath}/~{fname}.~{aligner}.bam
    rm ~{outputPath}/~{fname}.~{aligner}.sam
    fi
    
    if [ "~{aligner}" == "magic-blast" ];then
    zcat ~{outputPath}/~{fname}.R1.fq.gz | split -l4000000 -d --suffix-length 5 --filter='pigz -p 4 > $FILE.gz' - ~{outputPath}/~{fname}.R1.fq.split
    zcat ~{outputPath}/~{fname}.R2.fq.gz | split -l4000000 -d --suffix-length 5 --filter='pigz -p 4 > $FILE.gz' - ~{outputPath}/~{fname}.R2.fq.split
    declare -a filesR1=(~{outputPath}/~{fname}.R1.fq.split*)
    declare -a filesR2=(~{outputPath}/~{fname}.R2.fq.split*)

    do_magicblast()~{bracketOpen}
    magicblast \
    -db ~{index} \
    -query $1 \
    -query_mate $2 \
    -infmt fastq \
    -max_intron_length 500 \
    -outfmt sam > ~{outputPath}/~{fname}.split$3.sam
    ~{bracketClose}
    
    (
    for f in $~{bracketOpen}!filesR1[@]~{bracketClose};do
    ((i=i%~{threads}));((i++==0)) && wait
    do_magicblast "$~{bracketOpen}filesR1[$f]~{bracketClose}" "$~{bracketOpen}filesR2[$f]~{bracketClose}" "$f" &
    done && wait
    )
    
    (for f in $~{bracketOpen}!filesR1[@]~{bracketClose};do
    ((i=i%~{threads}));((i++==0)) && wait
    samtools view -b ~{outputPath}/~{fname}.split$f.sam | samtools sort -o ~{outputPath}/~{fname}.split$f.bam &
    done && wait
    )
    
    cmd="samtools merge -@ ~{threads} ~{outputPath}/~{fname}.magic-blast.bam"  
    for f in $~{bracketOpen}!filesR1[@]~{bracketClose}
    do
    cmd="$cmd ~{outputPath}/~{fname}.split$f.bam"
    done
    $cmd
    
    bash -c 'rm ~{outputPath}/~{fname}.split*'
    bash -c 'rm ~{outputPath}/~{fname}.R1.fq.split*'
    bash -c 'rm ~{outputPath}/~{fname}.R2.fq.split*'
    fi
    
  }
  
  output{
    File bamFile=outputPath+"/"+fname+"."+aligner+".bam"
  }
}

task filterBadAlignments{
  input{
    File bamFile
    String fname
    String outputPath
    String aligner
    Int threads
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    samtools view -@ ~{threads} -b -F 260 ~{bamFile} -o ~{outputPath}/~{fname}.~{aligner}.filter1.bam
    
    samtools view -@ ~{threads} ~{outputPath}/~{fname}.~{aligner}.filter1.bam \
    | awk 'sqrt($9*$9)>119 && sqrt($9*$9)<501' \
    | awk '~{bracketOpen}print $1~{bracketClose}' \
    > ~{outputPath}/~{fname}.~{aligner}.filter1.mapped_proper_pair.txt
    
    java -jar /usr/lib/picard.jar FilterSamReads \
    I=~{outputPath}/~{fname}.~{aligner}.filter1.bam \
    O=~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    READ_LIST_FILE=~{outputPath}/~{fname}.~{aligner}.filter1.mapped_proper_pair.txt \
    FILTER=includeReadList
    
    samtools view -@ ~{threads} ~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    | awk '~{bracketOpen}read=$0;sub(/.*NM:i:/,X,$0);sub(/\t.*/,X,$0);if(int($0)>7)~{bracketOpen}print read~{bracketClose}~{bracketClose}' \
    | awk '~{bracketOpen}print $1~{bracketClose}' \
    > ~{outputPath}/~{fname}.~{aligner}.filter2.high_mismatch.txt
    
    java -jar /usr/lib/picard.jar FilterSamReads \
    I=~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    O=~{outputPath}/~{fname}.~{aligner}.filter3.bam \
    READ_LIST_FILe=~{outputPath}/~{fname}.~{aligner}.filter2.high_mismatch.txt \
    FILTER=excludeReadList
    
  }
  
  output{
    File bamFilter=outputPath+"/"+fname+"."+aligner+".filter3.bam"
  }
}

task removeDuplicates{
  input{
    Boolean useUMI
    String bamFilter
    String fname
    String outputPath
    String aligner
    Int threads
  }
  
  #index bamFilter required by umi_tools; index filtered.dedup.bam for alignment visualization
  
  command{   
    if [[ ~{useUMI} == true ]];then
    samtools index -@ ~{threads} ~{bamFilter}
    
    umi_tools dedup --paired \
    -I ~{bamFilter} \
    -S ~{outputPath}/~{fname}.~{aligner}.filtered.dedup.bam \
    --output-stats=~{outputPath}/metrics/removeDuplicates-UMI_tools 
    
    else
    java -jar /usr/lib/picard.jar MarkDuplicates \
    I=~{bamFilter} \
    O=~{outputPath}/~{fname}.~{aligner}.filtered.dedup.bam \
    M=~{outputPath}/metrics/removeDuplicates-statsPicard.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=true
    
    fi
    
    samtools index -@ ~{threads} ~{outputPath}/~{fname}.~{aligner}.filtered.dedup.bam
  }
  
  output{
    File bamFilterDedup=outputPath+"/"+fname+"."+aligner+".filtered.dedup.bam"
  }
}


task getBamMetrics{
  input{
    String fname
    String outputPath
    String aligner
    String bamFilterDedup
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
    I=~{bamFilterDedup} \
    O=~{picardOut}/picard \
    VALIDATION_STRINGENCY=SILENT    
  }
  
  output{
    File picardInsertSizeMetrics=picardOut+'/picard.insert_size_metrics'
    File picardGcBiasMetrics=picardOut+'/picard.gc_bias.summary_metrics'
  }
  
}

task parseMethControl{
  input{
    String fname
    String outputPath
    String bamFilterDedup
    String aligner
    String? seqMeth
    String? seqUmeth
    Int threads
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    samtools view -@ ~{threads} ~{bamFilterDedup} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" ~{bracketOpen}print $2,$1~{bracketClose}' | sort -n -k1,1 > ~{outputPath}/metrics/meth_ctrl.counts
    total=$(samtools view -@ ~{threads} ~{bamFilterDedup} | wc -l)
    unmap=$(cat ~{outputPath}/meth_ctrl.counts | grep '^\*' | cut -f2); if [[ -z $unmap ]]; then unmap="0"; fi
    methyl=$(cat ~{outputPath}/meth_ctrl.counts | grep ~{seqMeth} | cut -f2); if [[ -z $methyl ]]; then methyl="0"; fi
    unmeth=$(cat ~{outputPath}/meth_ctrl.counts | grep ~{seqUmeth} | cut -f2); if [[ -z $unmeth ]]; then unmeth="0"; fi
    pct_meth_ctrl=$(echo "scale=3; ($methyl + $unmeth)/$total * 100" | bc -l); if [[ -z $pct_meth_ctrl ]]; then pct_meth_ctrl="0"; fi
    bet_meth_ctrl=$(echo "scale=3; $methyl/($methyl + $unmeth)" | bc -l); if [[ -z $bet_meth_ctrl ]]; then bet_meth_ctrl="0"; fi
    echo -e "total\tunmap\tmethyl\tunmeth\tPCT_METH_CTRL\tMETH_CTRL_BETA" > ~{outputPath}/metrics/meth_ctrl_summary.txt
    echo -e "$total\t$unmap\t$methyl\t$unmeth\t$pct_meth_ctrl\t$bet_meth_ctrl" >> ~{outputPath}/metrics/meth_ctrl_summary.txt
  }
  
  output{
    File methCtrlSummary=outputPath+"/metrics/meth_ctrl_summary.txt"
  }
}

task runMedips{
  input{
    File bamFilterDedup
    String fname
    String outputPath
    String aligner
    Int windowSize
  }
  
  String outMedips=outputPath+'/runMedips'
  
  command{
    mkdir -p ~{outMedips}
    
    r /workflow/runMedips/runMedips.r \
    --bamFile ~{bamFilterDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
    
    r /workflow/runMedips/runMedestrand.r \
    --bamFile ~{bamFilterDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
  }
  output{
    File medipsCount=outMedips+'/MEDIPS_hg38_'+fname+'_ws'+windowSize+'_count.txt.gz'
    File medipsRms=outMedips+'/MEDIPS_hg38_'+fname+'_ws'+windowSize+'_rms.txt.gz'
    File medestrandWig=outMedips+'/MeDESTrand_hg38_'+fname+'_ws'+windowSize+'_wig.bed.gz'
  }
}


