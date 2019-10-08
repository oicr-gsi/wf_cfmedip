# wf-cfMeDIP
Workflow for cfMeDIP data analysis using Docker

## User configuration: IMPORTANT!!! READ ME FIRST!!!

**Do not use sudo to execute `docker` commands**, instead, add an existing user to the _docker_ group `sudo usermod -aG docker harrycallahan`, which grants this user permissions to execute the `docker` command (i.e. `docker image ls`, `docker run`, `docker build`, etc).

To execute `docker run [docker_image]` in a development server, the user must pass its user and group IDs (in numeric format) to the `docker` command: `docker run --rm -u $(id -u):$(id -g) -ti [docker_image]`. As a consequence, the container will trigger processes in the host machine that belong to the user, making things like writing data to network shares possible. This user does not exist in the container, and the command prompt shows `I have no name!`, but processes are triggered under the valid UID nonetheless.

## Install & run
Download repository (currently private repo, Github user and password will be requested):
`git clone https://github.com/translational-genomics-laboratory/wf-cfMeDIP.git`

Build the image (can take more than an hour, requires sudo privileges in the host system): 
`sudo docker build -t wf_cfmedip:latest wf-cfMeDIP/`

Launch docker container in bash mode:
`docker run --rm -u $(id -u):$(id -g) -ti -v /host_folders/docker/:/home/docker/ wf_cfmedip:latest /bin/bash`

## Running jobs
The workflow is trigered by the script `wf_main.R` as follows:
`docker run --rm -u $(id -u):$(id -g) -v /host_folders/docker/:/home/docker/ wf_cfmedip:latest Rscript /home/R/wf_main.R [params]`.


