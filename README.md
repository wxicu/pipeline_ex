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

Dockerfile: 
  1. use Dockerfile to create a docker image containing demuxlet
  2. the image now is also available in the docker hub

sing_demuxlet.sif: 
  1. singularity image containing demuxlet (upload failure, try again later, maybe a little bit large)
  2. can use dockerfile from the docker hub directly
