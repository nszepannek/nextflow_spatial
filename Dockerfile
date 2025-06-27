FROM ubuntu:24.04

# Essentials
RUN apt-get update && apt-get install -y \
    curl \
    default-jre \
    r-base \
    git \
    wget \
    bzip2 \
    tar \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Download and extract Space Ranger v3.1.3
RUN curl -o spaceranger-3.1.3.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.3.tar.gz?Expires=1751018721&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=UjFUE9FZFi7hMN4drlDI3zOpfO9xnsPwW6DV1QpkDahRBf5yAd1LJLong6Lvf2GFFyjgd6vTe6EJ-Yy-KgcTPlo-7i-FKmcdGRiR8m~kTPUW2lqdTonJhUUERLBtnki-n7B~~YRM4k9HwyebBjajO7cDZvOi6N23fdFfrdPUUC9nG7ejzKH2y8W2UO1meiDRhLrhQa6yLhfwmDomJnXaLdy2Fy9ctv1siGHt4wM2f~iM7n0Led3uNyRHRysKLjSHly4ymZLPXUOWZu7SR4BpOHLbvR~MdIJ23LFjQr~Dy5EXIc6MjMNLHLQ62BRrD8cC7csQ3aUbq3KC5qKG0mtVqQ__" && \
    tar -xzf spaceranger-3.1.3.tar.gz && \
    rm spaceranger-3.1.3.tar.gz

# Add spaceranger to PATH
ENV PATH="/opt/spaceranger-3.1.3:$PATH"

# Install nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# set workdir
WORKDIR /workspace

# Kopiere dein Skript
COPY . /workspace/

# Installiere Pakete (einmalig beim Bauen)
RUN Rscript /workspace/scripts/install_and_run.R
