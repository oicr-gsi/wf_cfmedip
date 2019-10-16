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
		bowtie2 \
		samtools \
		trimmomatic

RUN curl -L -o /usr/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/2.20.8/picard.jar

RUN cd /home \
	&& wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2019-09-12.tar.gz \
	&& tar xvzf gmap-gsnap-2019-09-12.tar.gz \
	&& cd gmap-2019-09-12 \
	&& ./configure && make && make install \
	&& cd / \
	&& rm -rf /home/gmap-2019-09-12 && rm /home/gmap-gsnap-2019-09-12.tar.gz
	

RUN R -e 'install.packages(c("BiocManager","optparse","reshape2"))' \
	&& R -e 'library(BiocManager);BiocManager::install(c("MEDIPS","BSgenome.Hsapiens.UCSC.hg19"))'


RUN mkdir /home/data \
	&& mkdir /home/R

COPY data/*.gz /home/data/
COPY R/*.R /home/R/


ENTRYPOINT ["Rscript /home/R/wf_main.R"]
CMD [""]


