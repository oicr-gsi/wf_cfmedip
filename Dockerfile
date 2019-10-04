FROM r-base:3.6.1
    
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		libxml2-dev \
		libcurl4-openssl-dev \
		libssl-dev \
		default-jre \
		bowtie2 \
		samtools \
		trimmomatic


RUN R -e 'install.packages(c("BiocManager","optparse","reshape2"))' \
	R -e 'library(BiocManager);BiocManager::install(c("MEDIPS","BSgenome.Hsapiens.UCSC.hg19"))'


RUN mkdir /home/data \
	&& mkdir /home/R

COPY data/*.gz /home/data/
COPY R/*.R /home/R/


ENTRYPOINT ["Rscript /home/R/wf_main.R"]
CMD [""]

