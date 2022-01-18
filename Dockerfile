
### Some code adapted from 
### https://github.com/tschaffter/rstudio/blob/main/Dockerfile

# Base image
FROM rocker/tidyverse:latest

ENV miniconda3_version="py39_4.9.2"
ENV miniconda_bin_dir="/opt/miniconda/bin"
ENV PATH="${PATH}:${miniconda_bin_dir}"

RUN apt-get update -qq -y \
    && apt-get install --no-install-recommends -qq -y \
        bash-completion \
        curl \
        gosu \
        libxml2-dev \
        zlib1g-dev \
        libxtst6 \
        libxt6 \
        git \
    && apt-get -y autoclean \
    && apt-get -y autoremove \
    && curl -fsSLO https://repo.anaconda.com/miniconda/Miniconda3-${miniconda3_version}-Linux-x86_64.sh \
    && bash Miniconda3-${miniconda3_version}-Linux-x86_64.sh \
        -b \
        -p /opt/miniconda \
    && rm -f Miniconda3-${miniconda3_version}-Linux-x86_64.sh \
    && useradd -u 1500 -s /bin/bash miniconda \
    && chown -R miniconda:miniconda /opt/miniconda \
    && chmod -R go-w /opt/miniconda \
    && conda --version



# install R packages

# Install R packages
RUN install2.r --error \
    ggrepel \
    devtools \
    cowplot \
    VGAM \
    mclust \
    fpc \
    matrixStats \
    readr

# install R packages
RUN R -q -e "install.packages('BiocManager')"
RUN R -q -e "BiocManager::install(c('optparse','KernSmooth','ks','lattice','gridExtra'))"
RUN R -e "devtools::install_github('caravagnalab/BMix')"
RUN R -e "devtools::install_github('keyuan/ccube@6ad86869d8967ed4b1df25f4a9837c879c5b2601')"
RUN R -e "devtools::install_github('Wedge-Oxford/dpclust@75f5d7ef1e3e53585f86801fde76dd4c4aa86324')"


WORKDIR /letter


# Create conda environments

# create pyclone env
RUN conda init bash \
    && git clone https://github.com/Roth-Lab/pyclone-vi.git \
    && cd pyclone-vi \
    && conda create -c conda-forge -n pyclone-vi --file requirements.txt --yes \
    && rm -rf ./pyclone-vi


SHELL ["conda", "run", "-n", "pyclone-vi", "/bin/bash", "-c"]


RUN pip install git+https://github.com/Roth-Lab/pyclone-vi.git 

    





