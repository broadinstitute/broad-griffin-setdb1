############################################################
# Dockerfile for TE pipeline
# SmartMap
# Based on Debian slim
############################################################


FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

LABEL maintainer = "Eugenio Mattei"
LABEL software = "TE pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Repeats-quantification-chromatin"

ENV BOWTIE2_VERSION 2.4.3

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    build-essential \
    cpanminus \
    git \
    liblz4-dev \
    liblzma-dev \
    libncurses5-dev \
    libbz2-dev \
    unzip \
    wget \
    zlib1g-dev &&\
    rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
WORKDIR /software
ENV PATH="/software:${PATH}"

RUN cpanm Sys::Hostname

# Install Bowtie2 2.3.4.3
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-source.zip && \
    unzip bowtie2-${BOWTIE2_VERSION}-source.zip && cd bowtie2-${BOWTIE2_VERSION} && make static-libs && make STATIC_BUILD=1 && \
    cp bowtie2* .. && \
    cd .. && rm -rf bowtie2-${BOWTIE2_VERSION}*

RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1.zip && \
    unzip hisat2-2.2.1.zip && cd hisat2-2.2.1 && cp hisat2* .. && cd .. && \
    rm -rf hisat2-2.2.1

RUN wget https://github.com/shah-rohan/SmartMap/archive/refs/heads/master.zip && unzip master.zip && cd SmartMap-master && \
    make && mv SmartMap* .. && mv Default/SmartMap ../ && cd .. && rm -rf SmartMap-master


FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "Share-seq pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="SmartMap"

RUN apt-get update && apt-get install -y \
    gzip \
    cpanminus &&\
    rm -rf /var/lib/apt/lists/*

# Create and setup new user
ENV USER=shareseq
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"
RUN cpanm Sys::Hostname

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /software/bowtie2* /software/
COPY --from=builder --chown=$USER:$USER /software/hisat2* /software/
COPY --from=builder --chown=$USER:$USER /software/SmartMap* /software/
COPY --from=builder /usr/lib/x86_64-linux-gnu/perl/5.36 /usr/lib/x86_64-linux-gnu/perl/5.36/
COPY --from=builder --chown=$USER:$USER /lib/x86_64-linux-gnu/* /lib/x86_64-linux-gnu/
COPY --from=builder --chown=$USER:$USER /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6
#COPY --from=builder --chown=$USER:$USER /usr/local/bin/* /usr/local/bin/

#COPY --chown=$USER:$USER src/bash/monitor_script.sh /usr/local/bin

USER $USER
