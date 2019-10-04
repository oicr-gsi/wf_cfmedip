# wf-cfMeDIP
Workflow for cfMeDIP data analysis using Docker



## Install & run
Download repository:
`git clone https://github.com/translational-genomics-laboratory/wf-cfMeDIP.git`

Build the image (can take more than an hour); requires sudo privileges in the host system. 
`sudo docker build -t wf_cfmedip:latest wf-cfMeDIP/`

Launch docker container:
`docker run --rm -u $(id -u):$(id -g) -ti -v /data/results/cfMeDIP/docker/:/home/docker/ wf_cfmedip:latest`
The user inside the container should be a valid user in the host system (**not sudo or root user**); this makes the processes within the container being executed by a user with certain privileges in the host machine, and makes writing data to network shares possible.
