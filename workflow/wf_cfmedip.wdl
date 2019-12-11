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
      patternUMI2=patternUMI2
  }
  
  call alignReads {
    input:
      extrR1=extractUMI.extrR1,
      extrR2=extractUMI.extrR2,
      outputPath=outputPath,
      fname=fname,
      index=index,
      aligner=aligner
  }
  
  call filterBadAlignments{
    input:
      bamFile=alignReads.bamFile,
      outputPath=outputPath,
      fname=fname,
      aligner=aligner
  }
  
  call removeDuplicates{
    input:
      useUMI=useUMI,
      bamFilter=filterBadAlignments.bamFilter,
      outputPath=outputPath,
      fname=fname,
      aligner=aligner
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
      seqUmeth=seqUmeth
  }
  
  call getFilterMetrics{
    input:
      extrR1=extractUMI.extrR1,
      bamFilterDedup=removeDuplicates.bamFilterDedup,
      outputPath=outputPath,
      aligner=aligner,
      fname=fname
  }
  
  call doPicardDedup{
    input:
      bamFilter=filterBadAlignments.bamFilter,
      fname=fname,
      outputPath=outputPath,
      aligner=aligner
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
  }
  
  command{
    mkdir -p $(dirname ~{outputPath})
    
    if [[ ~{useUMI} == true ]];then  
    umi_tools extract --extract-method=string \
    --bc-pattern=~{patternUMI} --bc-pattern2=~{patternUMI2} -L extract.log \
    -I ~{R1} \
    --read2-in=~{R2} \
    -S ~{outputPath}/~{fname}.R1.fq.gz \
    --read2-out=~{outputPath}/~{fname}.R2.fq.gz
    else
      cp ~{R1} ~{outputPath}/~{fname}.R1.fq.gz
    cp ~{R2} ~{outputPath}/~{fname}.R2.fq.gz
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
  }
  
  command{
    if [ "~{aligner}" == "bowtie2" ];then
    bowtie2 -p 8 -x ~{index} \
    -1 ~{extrR1} \
    -2 ~{extrR2} \
    -S ~{outputPath}/~{fname}.bowtie2.sam
    fi
    
    if [ "~{aligner}" == "gsnap" ];then
    gsnap -t 8 -m 0.10 -A sam \
    --gunzip --trim-mismatch-score=0 \
    -D $(dirname ~{index}) \
    -d $(basename ~{index}) \
    ~{extrR1} \
    ~{extrR2} \
    > ~{outputPath}/~{fname}.gsnap.sam
    fi
    
    if [ "~{aligner}" == "bwa" ];then
    bwa mem -t 8 \
    ~{index} \
    ~{extrR1} \
    ~{extrR2} \
    > ~{outputPath}/~{fname}.bwa.sam
    fi
    
    samtools view -b ~{outputPath}/~{fname}.~{aligner}.sam | samtools sort -o ~{outputPath}/~{fname}.~{aligner}.bam
    
    rm ~{outputPath}/~{fname}.~{aligner}.sam
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
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    samtools view -b -F 260 ~{bamFile} -o ~{outputPath}/~{fname}.~{aligner}.filter1.bam
    
    samtools view ~{outputPath}/~{fname}.~{aligner}.filter1.bam \
    | awk 'sqrt($9*$9)>119 && sqrt($9*$9)<501' \
    | awk '~{bracketOpen}print $1~{bracketClose}' \
    > ~{outputPath}/~{fname}.~{aligner}.filter1.mapped_proper_pair.txt
    
    java -jar /usr/lib/picard.jar FilterSamReads \
    I=~{outputPath}/~{fname}.~{aligner}.filter1.bam \
    O=~{outputPath}/~{fname}.~{aligner}.filter2.bam \
    READ_LIST_FILE=~{outputPath}/~{fname}.~{aligner}.filter1.mapped_proper_pair.txt \
    FILTER=includeReadList
    
    samtools view ~{outputPath}/~{fname}.~{aligner}.filter2.bam \
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
  }
  
  command{
    if [[ ~{useUMI} == true ]];then
    samtools index ~{bamFilter}
    
    umi_tools dedup --paired \
    -I ~{bamFilter} \
    -S ~{outputPath}/~{fname}.~{aligner}.filtered.dedup.bam \
    --output-stats=deduplicated 
    else
      java -jar /usr/lib/picard.jar MarkDuplicates \
    I=~{bamFilter} \
    O=~{outputPath}/~{fname}.~{aligner}.filtered.dedup.bam \
    M=~{outputPath}/~{fname}.~{aligner}.filtered.dedup-statsPicard.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=true
    fi
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
    R=~{fasta} \
    I=~{bamFilterDedup}\
    O=~{picardOut}/~{fname}.~{aligner} \
    VALIDATION_STRINGENCY=SILENT
    
    java -jar /usr/lib/picard.jar CollectGcBiasMetrics \
    R=~{fasta} \
    I=~{bamFilterDedup} \
    O=~{picardOut}/~{fname}.~{aligner}.gc_bias_metrics.txt \
    S=~{picardOut}/~{fname}.~{aligner}.summary_gc_bias_metrics.txt \
    CHART=~{picardOut}/~{fname}.~{aligner}.gc_bias_metrics.pdf  
    
  }
  
  output{
    File picardMultipleMetrics=picardOut+'/'+fname+'.'+aligner+'.alignment_summary_metrics.txt'
    File picardGcBiasMetrics=picardOut+'/'+fname+'.'+aligner+'.gc_bias_metrics.txt'
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
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    samtools view ~{bamFilterDedup} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" ~{bracketOpen}print $2,$1~{bracketClose}' | sort -n -k1,1 > ~{outputPath}/meth_ctrl.counts
    total=$(samtools view ~{bamFilterDedup} | wc -l)
    unmap=$(cat ~{outputPath}/meth_ctrl.counts | grep '^\*' | cut -f2); if [[ -z $unmap ]]; then unmap="0"; fi
    methyl=$(cat ~{outputPath}/meth_ctrl.counts | grep ~{seqMeth} | cut -f2); if [[ -z $methyl ]]; then methyl="0"; fi
    unmeth=$(cat ~{outputPath}/meth_ctrl.counts | grep ~{seqUmeth} | cut -f2); if [[ -z $unmeth ]]; then unmeth="0"; fi
    pct_meth_ctrl=$(echo 'scale=3; ($methyl + $unmeth)/$total * 100' | bc -l); if [[ -z $pct_meth_ctrl ]]; then pct_meth_ctrl="0"; fi
    bet_meth_ctrl=$(echo 'scale=3; $methyl/($methyl + $unmeth)' | bc -l); if [[ -z $bet_meth_ctrl ]]; then bet_meth_ctrl="0"; fi
    echo -e "total\tunmap\tmethyl\tunmeth\tPCT_METH_CTRL\tMETH_CTRL_BETA" > ~{outputPath}/meth_ctrl_summary.txt
    echo -e "$total\t$unmap\t$methyl\t$unmeth\t$pct_meth_ctrl\t$bet_meth_ctrl" >> ~{outputPath}/meth_ctrl_summary.txt
  }
  
  output{
    File methCtrlSummary=outputPath+"/meth_ctrl_summary.txt"
  }
}

task runMedips{
  input{
    File bamFilterDedup
    String fname
    String outputPath
    String aligner
    String windowSize
  }
  
  String outMedips=outputPath+'/runMedips'
  
  command{
    mkdir -p ~{outMedips}
    
    Rscript /workflow/runMedips/runMedips.R \
    --bamFile ~{bamFilterDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
    
    Rscript /workflow/runMedips/runMedestrand.R \
    --bamFile ~{bamFilterDedup} \
    --outputDir ~{outMedips} \
    --windowSize ~{windowSize}
  }
  output{
    File medestrandWig=outMedips+'/MeDESTrand_hg38_'+fname+'_ws'+windowSize+'_wig.bed.gz'
  }
}

task getFilterMetrics{
  input{
    File extrR1
    File bamFilterDedup
    String outputPath
    String aligner
    String fname
  }
  
  command{
    total=`echo "$(gunzip -k -c ~{extrR1} | wc -l)/2" | bc`
    filter1=$(samtools view ~{outputPath}/~{fname}.~{aligner}.filter1.bam | wc -l)
    filter2=$(samtools view ~{outputPath}/~{fname}.~{aligner}.filter2.bam | wc -l)
    filter3=$(samtools view ~{outputPath}/~{fname}.~{aligner}.filter3.bam | wc -l)
    dedup=$(samtools view ~{bamFilterDedup} | wc -l)
    echo -e "total\tfilter1\tfilter2\tfilter3\tdedup" > ~{outputPath}/filter_metrics.txt
    echo -e "$total\t$filter1\t$filter2\t$filter3\t$dedup" >> ~{outputPath}/filter_metrics.txt
  }
  
  output{
    File filterMetrics=outputPath+"/filter_metrics.txt"
  }
}

task doPicardDedup{
  input{
    File bamFilter
    String fname
    String outputPath
    String aligner
  }
  
  String picardOut=outputPath+'/picard'
  
  command{
    mkdir -p ~{picardOut}
    
    java -jar /usr/lib/picard.jar MarkDuplicates \
    I=~{bamFilter} \
    O=~{outputPath}/~{fname}.~{aligner}.filter.dedup-Picard.bam \
    M=~{picardOut}/~{fname}.~{aligner}.filter.dedup-Picard.metrics \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=true
  }
}

