# pipeline_ex
pipeline exercise with nextflow 

main.nf:
  1. demuxlet

nextflow.config:
  1. parameter info
  2. docker
  3. singularity
  4. conda
  5. use profile parameter to switch between different containers

Dockerfile_old: 
  1. use Dockerfile to create a docker image containing old-version demuxlet 
 
Dockerfile_new: 
  1. use Dockerfile to create a docker image containing popscle as well as new-version demuxlet 
  2. available from the github 
  
Singularity:
  1. sif upload failure, try again later
  2. can also use docker image in docker hub directly
