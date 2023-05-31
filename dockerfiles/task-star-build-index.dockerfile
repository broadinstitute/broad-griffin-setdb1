############################################################
# Dockerfile for TE pipeline
# Based on Debian slim
############################################################

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f as builder

ENV STAR_VERSION 2.7.10b

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    build-essential \
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

# Compile STAR
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz && tar -xzf ${STAR_VERSION}.tar.gz && \cd STAR-${STAR_VERSION}/source && make STAR && rm ../../${STAR_VERSION}.tar.gz && mv /software/STAR-${STAR_VERSION}/bin/Linux_x86_64/* /usr/local/bin/

FROM debian@sha256:3ecce669b6be99312305bc3acc90f91232880c68b566f257ae66647e9414174f

LABEL maintainer = "Eugenio Mattei"
LABEL software = "TE pipeline"
LABEL software.version="0.0.1"
LABEL software.organization="Broad Institute of MIT and Harvard"
LABEL software.version.is-production="No"
LABEL software.task="STAR-index-build"

# Create and setup new user
ENV USER=te
WORKDIR /home/$USER

RUN groupadd -r $USER &&\
    useradd -r -g $USER --home /home/$USER -s /sbin/nologin -c "Docker image user" $USER &&\
    chown $USER:$USER /home/$USER

# Add folder with software to the path
ENV PATH="/software:${PATH}"

# Copy the compiled software from the builder
COPY --from=builder --chown=$USER:$USER /usr/local/bin/* /usr/local/bin/
COPY --from=builder --chown=$USER:$USER /usr/lib/x86_64-linux-gnu/libgomp.so.1 /lib/x86_64-linux-gnu/libncurses.so.6 /lib/x86_64-linux-gnu/

USER $USER
