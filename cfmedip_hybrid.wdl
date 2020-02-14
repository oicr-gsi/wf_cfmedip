version 1.0

import "tasks/extractUMI.wdl" as extractUMI
import "tasks/alignMagicblast.wdl" as alignMagicblast
import "tasks/filterAlignments.wdl" as filterAlignments
import "tasks/alignBowtie2.wdl" as alignBowtie2
import "tasks/coupleBam.wdl" as coupleBam
import "tasks/removeDuplicates.wdl" as removeDuplicates
import "tasks/parseMethControl.wdl" as parseMethControl
import "tasks/getBamMetrics.wdl" as getBamMetrics
import "tasks/runMedips.wdl" as runMedips

workflow cfmedip_hybrid {
  input{
    String outputPath
    String indexBowtie2
    String indexMagicblast
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
  
    call alignMagicblast.alignMagicblast as alignMagicblast {
    input:
      extrR1=extractUMI.extrR1,
      extrR2=extractUMI.extrR2,
      outputPath=outputPath,
      fname=fname,
      indexMagicblast=indexMagicblast,
      threads=threads
  }
  
    call filterAlignments.filterAlignments as filterAlignmentsMagicblast{
    input:
      bamAligned=alignMagicblast.bamAligned,
      outputPath=outputPath,
      fname=fname,
      aligner="magic-blast",
      maxMismatch=6,
      threads=threads
  }
  
    call alignBowtie2.alignBowtie2 as alignBowtie2 {
    input:
      extrR1=extractUMI.extrR1,
      extrR2=extractUMI.extrR2,
      outputPath=outputPath,
      fname=fname,
      indexBowtie2=indexBowtie2,
      threads=threads,
      processHold=alignMagicblast.bamAligned #holds bowtie2 until magicblast is done and threads are available
  }
  
  call filterAlignments.filterAlignments as filterAlignmentsBowtie2{
    input:
      bamAligned=alignBowtie2.bamAligned,
      outputPath=outputPath,
      fname=fname,
      aligner="bowtie2",
      maxMismatch=1,
      threads=threads
  }
  
  call coupleBam.coupleBam as coupleBam{
    input:
      bamFiltered1=filterAlignmentsMagicblast.bamFiltered,
      bamFiltered2=filterAlignmentsBowtie2.bamFiltered,
      outputPath=outputPath,
      fname=fname,
      threads=threads
  }
  
    call removeDuplicates.removeDuplicates as removeDuplicates{
    input:
      useUMI=useUMI,
      bamFiltered=coupleBam.bamHybrid,
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
      windowSize=windowSize
  }


}
