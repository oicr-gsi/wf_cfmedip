# wf-cfMeDIP
Workflow for cfMeDIP data analysis using Docker

![wf_cfmedip_overview](img/plot_wf_cfmedip_overview.png)

## User configuration: IMPORTANT!!! READ ME FIRST!!!

**Do not use sudo to execute `docker` commands**, instead, add an existing user to the _docker_ group `sudo usermod -aG docker harrycallahan`, which grants this user permissions to execute the `docker` command (i.e. `docker image ls`, `docker run`, `docker build`, etc).

To execute `docker run [docker_image]` in a development server, the user must pass its user and group IDs (in numeric format) by using the following `docker` command: `docker run --rm -u $(id -u):$(id -g) -ti [docker_image]`. As a consequence, the container processes belong to a user that is known by the host, allowing things like reading and writing data located in network shares. This user does not exist in the container, and the command prompt shows `I have no name!`, but processes are triggered under the valid UID nonetheless.

## Download Docker image
A pre-built image is hosted in Docker Hub (private repository, login required). First, execute `docker login`, and then `docker pull albertojleon/wf_cfmedip` to create a mirror image in the local registry.

## Build image
Download repository (currently private repo, Github user and password will be requested):
`git clone https://github.com/translational-genomics-laboratory/wf-cfMeDIP.git`

Build the image (can take more than an hour): 
`docker build -t wf_cfmedip:latest wf-cfMeDIP/`

## Basic container configuration
When triggering the workflow, it is necessary to bind the folders containing 1) the short reads, 2) pre-built genomic index, and 3) working directory:
`/home/docker/R/wf_main.R`, which gets triggered as the docker image spawns:
`docker run --rm -u $(id -u):$(id -g) -v /host_folders/docker/:/home/docker/ albertojleon/wf_cfmedip:latest [script params]`.

To access the 


