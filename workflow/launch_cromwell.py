#sourced from: wf_cfmedip_07_launch_cromwell.py

import subprocess
import os
import json
import argparse

str_desc="Launch stript for local execution of the workflow. It builds the inputs.json file and then executes 'cromwell run wf_cfmedip.wdl -i inputs.json'. No parameter completeness checks are preformed as this is better handled during WDL execution"
str_usage="python3 launch_cromwell.py [required_args] [optional args]"

parser = argparse.ArgumentParser(description=str_desc, usage=str_usage)
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')

required.add_argument("--R1", dest="R1", help="First mate Fastq file")
required.add_argument("--R2", dest="R2", help="Second mate Fastq file")
required.add_argument("--aligner", dest="aligner", help="bwa, bowtie2 or gsnap")
required.add_argument("--indexPath", dest="indexPath", help="Genome index, including methylation control sequences (when appropriate)")
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

args = parser.parse_args()

R1=args.R1
R2=args.R2
aligner=args.aligner
indexPath=args.indexPath
outputPath=args.outputPath
fastaFile=args.fastaFile

os.makedirs(outputPath, exist_ok=True)

wf_inputs = { 
  'wf_cfmedip.R1':R1,
  'wf_cfmedip.R2':R2,
  'wf_cfmedip.aligner':aligner,
  'wf_cfmedip.index':indexPath,
  'wf_cfmedip.fasta':fastaFile,
  'wf_cfmedip.outputPath':outputPath}

if bool(args.patternUMI):
  wf_inputs['wf_cfmedip.patternUMI']=args.patternUMI

if bool(args.patternUMI2):
  wf_inputs['wf_cfmedip.patternUMI2']=args.patternUMI2

if bool(args.seqMeth):
  wf_inputs['wf_cfmedip.seqMeth']=args.seqMeth

if bool(args.seqUmeth):
  wf_inputs['wf_cfmedip.seqUmeth']=args.seqUmeth

if bool(args.useUMI):
  wf_inputs['wf_cfmedip.useUMI']=args.useUMI

if bool(args.windowSize):
  wf_inputs['wf_cfmedip.windowSize']=args.windowSize

if bool(args.sampleName):
  wf_inputs['wf_cfmedip.sampleName']=args.sampleName


with open(outputPath+'/wf_cfmedip.inputs.json', 'w') as json_file:
  json.dump(wf_inputs,json_file,indent=4)


#-Duser.dir sets java's working directory, used by cromwell to store execution files (without this, cromwell may use the user's home directory)
cmd = ['java','-Xmx1g','-Duser.dir='+outputPath,'-jar','/usr/lib/cromwell.jar',
        'run','/workflow/wf_cfmedip.wdl',
        '-i',outputPath+'/wf_cfmedip.inputs.json'
      ]

print(cmd)

proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
o, e = proc.communicate()
print('Output: ' + o.decode('ascii'))
print('Error: '  + e.decode('ascii'))
print('code: ' + str(proc.returncode))


