FROM ubuntu
RUN apt-get update && apt-get install -y\
	git \
	build-essential \
	libtool \
	wget \
	autoconf \
	automake \
	make \
	gcc \
	zlib1g-dev \
	libbz2-dev \	
	liblzma-dev \
	libcurl4-gnutls-dev \
	libssl-dev \
	libncurses5-dev


RUN git clone https://github.com/samtools/htslib \ 
	&& cd htslib \ 
	&& git submodule update --init --recursive \
	&& autoreconf -i \
	&& ./configure \
	&& make \
	&& make install \
	&& cd ..

RUN git clone https://github.com/statgen/demuxlet.git \
	&& cd demuxlet \
	&& autoreconf -vfi \
	&& ./configure \
	&& make \
	&& make install


