FROM python:3.8-slim-buster

ARG USERNAME=ovation
ARG USER_UID=1015
ARG USER_GID=$USER_UID

# Add ovation user/group
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME

# Install PR0J for Ovation and wget for conda
RUN apt-get update \
    && apt-get -y install proj-bin wget gcc \

ENV PATH /opt/conda/bin:$PATH

# basemap can't be installed without conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Basemap
RUN conda install -y -c anaconda basemap
ENV PROJ_LIB=/usr/share/proj

RUN pip --disable-pip-version-check --no-cache-dir requests install matplotlib numpy scipy aacgmv2 image

# make necessary directories
RUN mkdir Output
RUN chown ovation:ovation Output

USER $USERNAME

# Add the source code into the image
COPY --chown=ovation:ovation source source

COPY --chown=ovation:ovation SW_Data SW_Data
COPY --chown=ovation:ovation configuration configuration
COPY --chown=ovation:ovation source_plots source_plots

WORKDIR source
