version 1.0

task alignMagicblast{
  input{
    String extrR1
    String extrR2
    String indexMagicblast
    String fname
    String outputPath
    Int threads
  }
  
  String aligner="magic-blast"
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    zcat ~{outputPath}/~{fname}.R1.fq.gz | split -l4000000 -d --suffix-length 5 --filter='pigz -p 4 > $FILE.gz' - ~{outputPath}/~{fname}.R1.fq.split
    zcat ~{outputPath}/~{fname}.R2.fq.gz | split -l4000000 -d --suffix-length 5 --filter='pigz -p 4 > $FILE.gz' - ~{outputPath}/~{fname}.R2.fq.split
    declare -a filesR1=(~{outputPath}/~{fname}.R1.fq.split*)
    declare -a filesR2=(~{outputPath}/~{fname}.R2.fq.split*)
    
    do_magicblast()~{bracketOpen}
    magicblast \
    -db ~{indexMagicblast} \
    -query $1 \
    -query_mate $2 \
    -infmt fastq \
    -max_intron_length 500 \
    -outfmt sam > ~{outputPath}/~{fname}.~{aligner}.split$3.sam
    ~{bracketClose}
    
    (
      for f in $~{bracketOpen}!filesR1[@]~{bracketClose};do
      ((i=i%~{threads}));((i++==0)) && wait
      do_magicblast "$~{bracketOpen}filesR1[$f]~{bracketClose}" "$~{bracketOpen}filesR2[$f]~{bracketClose}" "$f" &
        done && wait
    )
    
    (for f in $~{bracketOpen}!filesR1[@]~{bracketClose};do
      ((i=i%~{threads}));((i++==0)) && wait
      samtools view -b ~{outputPath}/~{fname}.~{aligner}.split$f.sam | samtools sort -o ~{outputPath}/~{fname}.~{aligner}.split$f.bam &
        done && wait
    )
    
    cmd="samtools merge -@ ~{threads} ~{outputPath}/~{fname}.~{aligner}.bam"  
    for f in $~{bracketOpen}!filesR1[@]~{bracketClose}
    do
    cmd="$cmd ~{outputPath}/~{fname}.~{aligner}.split$f.bam"
    done
    $cmd
    
    bash -c 'rm ~{outputPath}/~{fname}.~{aligner}.split*'
    bash -c 'rm ~{outputPath}/~{fname}.R1.fq.split*'
    bash -c 'rm ~{outputPath}/~{fname}.R2.fq.split*'
    
  }
  
  output{
    File bamAligned=outputPath+"/"+fname+"."+aligner+".bam"
  }
}