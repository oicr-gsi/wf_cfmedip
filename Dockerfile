#rocker/r-ver is built on debian:stable, whereas r-base follows debian:testing
FROM rocker/r-ver:3.6.2 

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		bc \
		bowtie2 \
		bzip2 \
		curl \
		default-jre \
		fastp \
		libbz2-dev \
		libcurl4-openssl-dev \
		liblzma-dev \
		libncurses5-dev \
		libssl-dev \
		libxml2-dev \
		nano \
		pigz \
		wget \
		zlib1g-dev

RUN curl -L -o /usr/lib/picard.jar https://github.com/broadinstitute/picard/releases/download/2.21.8/picard.jar \
	&& curl -L -o /usr/lib/cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/48/cromwell-48.jar

RUN cd /home \
	&& wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
	&& tar -xvjf bwa-0.7.17.tar.bz2 \
	&& cd bwa-0.7.17 \
	&& make \
	&& cd /home && mv /home/bwa-0.7.17 /usr/lib/ && rm bwa-0.7.17.tar.bz2
	
RUN cd /home \
	&& wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/magicblast/1.5.0/ncbi-magicblast-1.5.0-x64-linux.tar.gz \
	&& tar -xvzf ncbi-magicblast-1.5.0-x64-linux.tar.gz \
	&& cp ncbi-magicblast-1.5.0/bin/magicblast /usr/bin/ \
	&& cp ncbi-magicblast-1.5.0/bin/makeblastdb /usr/bin/ \
	&& rm -rf ncbi-magicblast-1.5.0-x64-linux.tar.gz

RUN cd /home \
	&& wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
	&& tar -xvjf samtools-1.10.tar.bz2 \
	&& cd samtools-1.10 \
	&& make && make install \
	&& cd .. \
	&& rm -rf samtools-1.10 && rm samtools-1.10.tar.bz2

RUN R -e 'install.packages(c("BiocManager","docopt","reshape2","remotes"))' \
	&& R -e 'library(BiocManager);BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38","MEDIPS"))' \
	&& R -e 'library(remotes);remotes::install_github("oicr-gsi/modelTsne")' \
	&& R -e 'library(remotes);remotes::install_github("jxu1234/MeDEStrand")'

RUN apt-get install -y r-cran-littler

#Flag '--no-install-recommends' unsuitable for python as it skips installation of required libraries
#Required before UMI-tools(must be present in advance): pip3 install Cython
RUN apt-get install -y python3-pip \
	&& pip3 install Cython \
	&& pip3 install UMI-tools 

#Does it prevent loss of $PATH in cluster?
RUN echo "export PATH=$PATH:/usr/bin:/usr/local/bin:/usr/lib/bwa-0.7.17" >> /etc/bash.bashrc \
	&& cat /etc/profile | awk '{ if ($0=="export PATH") print "export PATH=$PATH:/usr/lib/bwa-0.7.17"; else print $0}' > /etc/profile_tmp \
	&& mv /etc/profile_tmp /etc/profile

COPY workflow /workflow

ENV PATH=/usr/lib/bwa-0.7.17:$PATH
