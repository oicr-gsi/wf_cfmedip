# wf-cfMeDIP
Workflow for cfMeDIP data analysis using Docker

## User configuration: IMPORTANT!!! READ ME FIRST!!!
Having sudo privileges in the host system is required to build the image. During this process, it is the conatiner's root user the one executing the Dockerfile commands (i.e. apt-get, configure environment variables, etc).

**Do not use sudo to execute `docker` commands**, instead, add an existing user to the _docker_ group `sudo usermod -aG docker harrycallahan`, which grants this user permissions to execute the `docker` command (i.e. `harrycallahan@devmachine:/$ docker image ls`).

To properly run a container, the user must pass its user and group IDs (in numeric format) to the `docker` command: `docker run --rm -u $(id -u):$(id -g) -ti r-base`. As a consequence, the container will trigger processes in the host machine that belong to the user, making things like writing data to network shares possible. This user does not exist in the container, and the command prompt shows `I have no name!`, but processes are triggered under the valid UID nonetheless.

## Install & run
Download repository (currently private repo, Github user and password will be requested):
`git clone https://github.com/translational-genomics-laboratory/wf-cfMeDIP.git`

Build the image (can take more than an hour, requires sudo privileges in the host system): 
`sudo docker build -t wf_cfmedip:latest wf-cfMeDIP/`

Launch docker container:
`docker run --rm -u $(id -u):$(id -g) -ti -v /host_folders/docker/:/home/docker/ wf_cfmedip:latest /bin/bash`


