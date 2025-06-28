FROM ubuntu:24.04

# Essentials
RUN apt-get update && apt-get install -y \
    curl \
    graphviz \
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
    libhdf5-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt

# Download and extract Space Ranger v3.1.3
# Update this link on docker re-build https://www.10xgenomics.com/support/software/space-ranger/downloads/previous-versions
RUN curl -o spaceranger-3.1.3.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.3.tar.gz?Expires=1751144007&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=btn1D62ZaN5G6dphUlnIHFwAT4HhZKsZ8xuX1PwWWf~MmeAg8pwH~NZ78m8winOpbmhTpTSCAHtvQOJaRmwxD~pK8GzV~0Zg2KAiodgmt8kTHfFO8kCdAStGwFkYt-ZBtTVjMhE5kUrcT3lkSrpfNB8tlka9M5p4PDeR8bG-UN-4J8pu0RJKecj4Ul-w7ts6qRVG~YiMRbFmmZhHEu2KmhvMtDqg0XMRQ-i9XQtnjKcHHlMXAp7glF7urP~Iq98ml6i31QDHpmtRF~gA1U564JPMkxy2FpcbKlGA0eOK0986SfAXO~E0RNRyLxdKjVv3nHuccpO5XhPKw2Wb~lM19g__" && \
    tar -xzf spaceranger-3.1.3.tar.gz && \
    rm spaceranger-3.1.3.tar.gz

# Add spaceranger to PATH
ENV PATH="/opt/spaceranger-3.1.3:$PATH"

# Install nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

WORKDIR /workspace
COPY Packages.R /workspace/Packages.R
RUN Rscript /workspace/Packages.R
