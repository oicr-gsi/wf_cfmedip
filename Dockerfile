FROM r-base:3.6.1
    
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		curl \
		nano \
		gawk \
		libxml2-dev \
		libcurl4-openssl-dev \
		libssl-dev \
		default-jre \
		bowtie2

RUN curl -L -o /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/2.20.8/picard.jar
RUN curl -L -o /usr/bin/cromwell.jar https://github.com/broadinstitute/cromwell/releases/download/47/cromwell-47.jar
RUN curl -L -o /usr/bin/womtool.jar https://github.com/broadinstitute/cromwell/releases/download/47/womtool-47.jar

RUN cd /home \
	&& wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2019-09-12.tar.gz \
	&& tar xvzf gmap-gsnap-2019-09-12.tar.gz \
	&& cd gmap-2019-09-12 \
	&& ./configure && make && make install \
	&& cd .. \
	&& rm -rf /home/gmap-2019-09-12 && rm /home/gmap-gsnap-2019-09-12.tar.gz
	
RUN cd /home \
	&& wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
	&& tar -xvjf samtools-1.9.tar.bz2 \
	&& cd samtools-1.9 \
	&& make && make install \
	&& cd .. \
	&& rm -rf samtools-1.9 && rm samtools-1.9.tar.bz2
	
RUN R -e 'install.packages(c("BiocManager","optparse","reshape2,devtools"))' \
	&& R -e 'library(BiocManager);BiocManager::install(c("MEDIPS","BSgenome.Hsapiens.UCSC.hg38"))' \
	&& R -e 'devtools::install_github("jxu1234/MeDEStrand")'
	
RUN apt-get install -y --no-install-recommends python3-pip \
	&& pip3 install UMI-tools
	
RUN mkdir /home/data \
	&& mkdir /home/R

COPY data/*.gz /home/data/
COPY R/*.R /home/R/

ENTRYPOINT ["Rscript /home/R/wf_main.R"]
CMD [""]

