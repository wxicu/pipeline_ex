FROM ubuntu:latest

# apt update and install global requirements
RUN apt-get update && apt-get -y install freebayes

# apt clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*




