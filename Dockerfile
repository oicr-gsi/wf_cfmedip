#rocker/r-ver is built on debian:stable, whereas r-base follows debian:testing
FROM rocker/r-ver:3.6.1 
    
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		bc \
		bowtie2 \
		bzip2 \
		curl \
		default-jre \
		libbz2-dev \
		libcurl4-openssl-dev \
		liblzma-dev \
		libssl-dev \
		libxml2-dev \
		nano \
		wget \
		zlib1g-dev

RUN curl -L -o /usr/lib/picard.jar https://github.com/broadinstitute/picard/releases/download/2.20.8/picard.jar \
	&& curl -L -o /usr/lib/cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/47/cromwell-47.jar \
	&& curl -L -o /usr/lib/womtool.jar https://github.com/broadinstitute/cromwell/releases/download/47/womtool-47.jar

RUN cd /home \
	&& wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2019-09-12.tar.gz \
	&& tar xvzf gmap-gsnap-2019-09-12.tar.gz \
	&& cd gmap-2019-09-12 \
	&& ./configure && make && make install \
	&& cd /home && rm -rf /home/gmap-2019-09-12 && rm /home/gmap-gsnap-2019-09-12.tar.gz
	
RUN cd /home \
	&& wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
	&& tar -xvjf bwa-0.7.17.tar.bz2 \
	&& cd bwa-0.7.17 \
	&& make \
	&& cd /home && mv /home/bwa-0.7.17 /usr/lib/ && rm bwa-0.7.17.tar.bz2

RUN cd /home \
	&& wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
	&& tar -xvjf samtools-1.9.tar.bz2 \
	&& cd samtools-1.9 \
	&& make && make install \
	&& cd .. \
	&& rm -rf samtools-1.9 && rm samtools-1.9.tar.bz2
	
RUN R -e 'install.packages(c("BiocManager","optparse","reshape2","devtools"))' \
	&& R -e 'library(BiocManager);BiocManager::install(c("MEDIPS","BSgenome.Hsapiens.UCSC.hg38"))' \
	&& R -e 'library(devtools);devtools::install_github("jxu1234/MeDEStrand")'

#Flag '--no-install-recommends' unsuitable for python as it skips installation of required libraries
RUN apt-get install -y python3-pip \
	&& pip3 install UMI-tools 

#Does it prevent loss of $PATH in cluster?
RUN echo "export PATH=$PATH:/usr/bin:/usr/local/bin:/usr/lib/bwa-0.7.17" > /etc/profile.d/custom_paths.sh
	
COPY R /home
	

