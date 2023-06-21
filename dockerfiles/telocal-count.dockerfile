############################################################
# Dockerfile for TE pipeline
# Based on Debian slim
############################################################


FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "TE pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="Repeats-quantification"

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    unzip \
    wget &&\
    rm -rf /var/lib/apt/lists/*


RUN wget -q https://github.com/mhammell-laboratory/TEtranscripts/archive/master.zip && \
    unzip master.zip && \
    python3 -m pip install --break-system-packages pysam argparse && \
    cd TEtranscripts-master && \
    python3 setup.py install && \
    cd ../ && rm -rf TEtranscripts-master && rm master.zip

 RUN wget -q https://github.com/mhammell-laboratory/TElocal/archive/refs/heads/master.zip && \
     unzip master.zip && \
     cd TElocal-master && \
     python3 setup.py install && \
     cd ../ && rm -rf TElocal-master

# Create and setup new user
ENV USER=te
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

USER $USER

