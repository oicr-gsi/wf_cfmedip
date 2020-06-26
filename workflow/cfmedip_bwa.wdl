version 1.0

import "tasks/extractUMI.wdl" as extractUMI
import "tasks/alignBwa.wdl" as alignBwa
import "tasks/filterAlignments.wdl" as filterAlignments
import "tasks/removeDuplicates.wdl" as removeDuplicates
import "tasks/parseMethControl.wdl" as parseMethControl
import "tasks/getBamMetrics.wdl" as getBamMetrics
import "tasks/runMedips.wdl" as runMedips
import "tasks/runMedestrand.wdl" as runMedestrand

workflow cfmedip_bwa {
  input{
    String outputPath
    String indexBwa
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
    Boolean useMedestrand = false
    File ROIFile = "/data/UCSC-hg38-CpG.bed"
  }
  
  String fname=if defined(sampleName) then select_first([sampleName,""]) else sub(basename(R1),"(\.fq)?(\.fastq)?(\.gz)?", "")
  
  call extractUMI.extractUMI as extractUMI {
    input:
      useUMI=useUMI,
      R1=R1,
      R2=R2,
      outputPath=outputPath,
      fname=fname,
      patternUMI=patternUMI,
      patternUMI2=patternUMI2,
      newReadLen=newReadLen
  }
  
  call alignBwa.alignBwa as alignBwa {
    input:
      extrR1=extractUMI.extrR1,
      extrR2=extractUMI.extrR2,
      outputPath=outputPath,
      fname=fname,
      indexBwa=indexBwa,
      threads=threads
  }
  
  call filterAlignments.filterAlignments as filterAlignmentsBwa{
    input:
      bamAligned=alignBwa.bamAligned,
      outputPath=outputPath,
      fname=fname,
      aligner="bwa",
      maxMismatch=6,
      threads=threads
  }
  
  call removeDuplicates.removeDuplicates as removeDuplicates{
    input:
      useUMI=useUMI,
      bamFiltered=filterAlignmentsBwa.bamFiltered,
      outputPath=outputPath,
      fname=fname,
      threads=threads
  }
  
  call getBamMetrics.getBamMetrics as getBamMetrics{
    input:
      fname=fname,
      outputPath=outputPath,
      bamDedup=removeDuplicates.bamDedup,
      fasta=fasta
  }
  
  call parseMethControl.parseMethControl as parseMethControl{
    input:
      bamDedup=removeDuplicates.bamDedup,
      outputPath=outputPath,
      fname=fname,
      seqMeth=seqMeth,
      seqUmeth=seqUmeth,
      threads=threads
  }
  
  call runMedips.runMedips as runMedips{
    input:
      bamDedup=removeDuplicates.bamDedup,
      fname=fname,
      outputPath=outputPath,
      windowSize=windowSize,
      ROIFile=ROIFile
  }
  
  call runMedestrand.runMedestrand as runMedestrand{
    input:
      bamDedup=removeDuplicates.bamDedup,
      fname=fname,
      outputPath=outputPath,
      windowSize=windowSize,
      useMedestrand=useMedestrand
  }
  
  
}
