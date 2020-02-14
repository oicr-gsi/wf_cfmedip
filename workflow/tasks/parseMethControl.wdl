version 1.0

task parseMethControl{
  input{
    String fname
    String outputPath
    String bamDedup
    String? seqMeth
    String? seqUmeth
    Int threads
  }
  
  String bracketOpen="{"
  String bracketClose="}"
  
  command{
    samtools view -@ ~{threads} ~{bamDedup} | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" ~{bracketOpen}print $2,$1~{bracketClose}' | sort -n -k1,1 > ~{outputPath}/metrics/meth_ctrl.counts
    total=$(samtools view -@ ~{threads} ~{bamDedup} | wc -l)
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