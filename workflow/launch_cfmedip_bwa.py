import subprocess
import os
import json
import argparse

str_desc="Launch stript for local execution of the workflow. It builds the inputs.json file and then executes 'cromwell run cfmedip_bwa.wdl -i inputs.json'. No parameter completeness checks are preformed as this is better handled during WDL execution"
str_usage="python3 launch_cromwell.py [required_args] [optional args]"

parser = argparse.ArgumentParser(description=str_desc, usage=str_usage)
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')

required.add_argument("--R1", dest="R1", help="First mate Fastq file")
required.add_argument("--R2", dest="R2", help="Second mate Fastq file")
required.add_argument("--indexBwa", dest="indexBwa", help="Genome index built for Bwa, including methylation control sequences (when appropriate)")
required.add_argument("--fastaFile",dest="fastaFile",help="Reference genome sequence in fasta format, including methylation control sequences (when appropriate)")
required.add_argument("--outputPath", dest="outputPath", help="Destination folder")

optional.add_argument("--sampleName", dest="sampleName", help="If provided, used as basename, otherwhise, basename from --R1 is used")

#Default values are defined in the WDL file
optional.add_argument("--patternUMI", dest="patternUMI", help="umi_tools parameter --bc-pattern (default=NNNNN)")
optional.add_argument("--patternUMI2", dest="patternUMI2", help="umi_tools parameter --bc-pattern2 (default=NNNN)")
optional.add_argument("--seqMeth", dest="seqMeth", help="Name of methylated control sequence (default=F19K16)")
optional.add_argument("--seqUmeth", dest="seqUmeth", help="Name of unmethylated control sequence (default=F24B22)")
optional.add_argument("--useUMI", dest="useUMI", help="Do reads include UMI sequences? 'true' or 'false' (default=true)")
optional.add_argument("--windowSize", dest="windowSize", help="Genomic window size (bp) (default=200)")
optional.add_argument("--threads", dest="threads", help="Number of threads used by the aligner")
optional.add_argument("--newReadLen", dest="newReadLen", help="fastp --max_len1 and --max_len2. Default = -1 (trimming disabled)")
optional.add_argument("--useMedestrand", dest="useMedestrand", help="Perform MeDEStrand analysis? 'true' or 'false' (default = false)")
optional.add_argument("--ROIFile", dest="ROIFile", help="File with regions of interest for MeDIPs. Default = /data/UCSC-hg38-CpG.bed")

args = parser.parse_args()
outputPath=args.outputPath

os.makedirs(outputPath, exist_ok=True)

wf_inputs = { 
  'cfmedip_bwa.R1':args.R1,
  'cfmedip_bwa.R2':args.R2,
  'cfmedip_bwa.indexBwa':args.indexBwa,
  'cfmedip_bwa.fasta':args.fastaFile,
  'cfmedip_bwa.outputPath':outputPath}

if bool(args.patternUMI):
  wf_inputs['cfmedip_bwa.patternUMI']=args.patternUMI

if bool(args.patternUMI2):
  wf_inputs['cfmedip_bwa.patternUMI2']=args.patternUMI2

if bool(args.seqMeth):
  wf_inputs['cfmedip_bwa.seqMeth']=args.seqMeth

if bool(args.seqUmeth):
  wf_inputs['cfmedip_bwa.seqUmeth']=args.seqUmeth

if bool(args.useUMI):
  wf_inputs['cfmedip_bwa.useUMI']=args.useUMI

if bool(args.windowSize):
  wf_inputs['cfmedip_bwa.windowSize']=args.windowSize

if bool(args.sampleName):
  wf_inputs['cfmedip_bwa.sampleName']=args.sampleName

if bool(args.threads):
  wf_inputs['cfmedip_bwa.threads']=args.threads

if bool(args.newReadLen):
  wf_inputs['cfmedip_bwa.newReadLen']=args.newReadLen

if bool(args.useMedestrand):
  wf_inputs['cfmedip_bwa.useMedestrand']=args.useMedestrand

if bool(args.ROIFile):
  wf_inputs['cfmedip_bwa.useROIFile']=args.ROIFile

with open(outputPath+'/cfmedip_bwa.inputs.json', 'w') as json_file:
  json.dump(wf_inputs,json_file,indent=4)


#-Duser.dir sets java's working directory, used by cromwell to store execution files (without this, cromwell may use the user's home directory)
cmd = ['java','-Xmx1g','-Duser.dir=/cromwell','-jar','/usr/lib/cromwell.jar',
      'run',
      '-i',outputPath+'/cfmedip_bwa.inputs.json',
      '/TGL/gsi/data/iScan/cfMeDIP/workflow/cfmedip_bwa.wdl'
      #'run','/workflow/cfmedip_bwa.wdl'
      
      ]

print(cmd)

proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
o, e = proc.communicate()
print('Output: ' + o.decode('ascii'))
print('Error: '  + e.decode('ascii'))
print('code: ' + str(proc.returncode))
