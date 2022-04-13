# Define base image
From continuumio/miniconda3

# Set working directory for the project, where is the yml file
WORKDIR /home/xic/docker_conda
 
# Create Conda environment from the YML file
COPY my_env.yml .
RUN conda env create -f my_env.yml

#solution1:
# Override default shell and use bash on the new environment
# SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]
# The code to run when container is started:
#ENTRYPOINT ["conda", "run", "-n", "myenv", "popscle", "--help"]

#Solution2:
# Activate the environment
RUN echo "source activate my_env" >> ~/.bashrc
ENV PATH /opt/conda/envs/myenv/bin:$PATH
#check with docker run -it conda:latest conda info --envs

# Demonstrate the environment is activated:
RUN echo "Make sure popscle is installed:"
RUN echo "popscle --help" >> ~/.bashrc


