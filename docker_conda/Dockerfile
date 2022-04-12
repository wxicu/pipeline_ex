# Define base image
From continuumio/miniconda3

# Set working directory for the project
WORKDIR /home/xic/docker_conda
 
# Create Conda environment from the YML file
COPY my_env.yml .
RUN conda env create -f my_env.yml
 
# Override default shell and use bash on the new environment
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]
 
# Demonstrate the environment is activated:
RUN echo "Make sure popscle is installed:"
RUN popscle --help

# The code to run when container is started:
#ENTRYPOINT ["conda", "run", "-n", "myenv", "docker", "pull", "ubuntu"]

