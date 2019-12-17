# wf-cfMeDIP
Workflow for cfMeDIP data analysis using Docker

![wf_cfmedip_overview](img/plot_wf_cfmedip_overview.png)

## User configuration: IMPORTANT!!! READ ME FIRST!!!

**Do not use sudo to execute `docker` commands**, instead, add an existing user to the _docker_ group `sudo usermod -aG docker harrycallahan`, which grants this user permissions to execute the `docker` command (i.e. `docker image ls`, `docker run`, `docker build`, etc).

To execute `docker run [docker_image]` in a development server, the user must pass its user and group IDs (in numeric format) by using the following `docker` command: `docker run --rm -u $(id -u):$(id -g) -ti [docker_image]`. As a consequence, the container processes belong to a user that is known by the host, allowing things like reading and writing data located in network shares. This user does not exist in the container, and the command prompt shows `I have no name!`, but processes are triggered under the valid UID nonetheless.

## Download pre-built image
A pre-built Docker image is hosted in Docker Hub. Execute `docker pull oicrgsi/wf_cfmedip` to create a mirror image in the local registry.

## Build image locally
Download repository:
`git clone https://github.com/oicr-gsi/wf-cfMeDIP.git`

Build the image (this process can take more than two hours): 
`docker build -t wf_cfmedip:latest wf-cfMeDIP/`

## Workflow parameters
| Parrameter  | Required/Optional | Description |
| --- | --- | --- |
| --R1  | Required | fastq file mate 1 (.gz allowed) |
| --R2  | Required | fastq file mate 2 (.gz allowed) |
| --aligner | Required | bowtie2 or bwa or gsnap |
| --index | Required | genomic index pre-built for the selected aligner |
| --fasta | Required | reference genomic sequence in fasta format |
| --outputPath | Required | output folder |
| --sampleName | Optional | labels files. If not provided, the filename from --R1 is used |
| --patternUMI | Optional | UMI-tools patternUMI parameter. Default "NNNNN" |
| --patternUMI2 | Optional | UMI-tools patternUMI2 parameter. Default = "NNNNN" |
| --seqMeth | Optional | Name of sequence to be used as methylated control. Default = "F19K16" |
| --seqUmeth | Optional | Name of sequence to be used as unmethylated control. Default = "F24B22" |
| --useUMI | Optional | Do reads contain UMIs? Default = true |
| --windowSize | Optional | MeDIPs window size parameter. Default = 200 |

## Example
![wf_cfmedip_file_structure](img/screenshot_file_structure.png)

```
|storage
	|cromwell
	|MyStudy
		|fastq
			|sample001_R1.fastq.gz
			|sample001_R2.fastq.gz
		|output
			|sample001
				|run_inside_container.sh
		|index
			|bowtie2
				|ucsc.hg38_F19K16_F24B22
					|ucsc.hg38_F19K16_F24B22
			|fasta
				|ucsc.hg38_F19K16_F24B22.fasta
```


The command:
```
python3 /workflow/launch_cromwell.py \
	--R1 /data/fastq/sample001_R1_001.fastq.gz \
	--R2 /data/fastq/sample001_R2_001.fastq.gz \
	--aligner bowtie2 \
	--indexPath /data/index/bowtie2/ucsc.hg38_F19K16_F24B22/ucsc.hg38_F19K16_F24B22 \
	--fastaFile /data/index/fasta/ucsc.hg38_F19K16_F24B22.fasta \
	--outputPath /data/output/sample001
```









