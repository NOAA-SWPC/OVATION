FROM continuumio/anaconda3

ARG USERNAME=ovation
ARG USER_UID=1015
ARG USER_GID=$USER_UID

# Add ovation user/group
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME

# Install PR0J for Ovation
RUN apt-get update \
    && apt-get -y install proj-bin gcc \
# Clean up
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Basemap
RUN conda install -y -c anaconda basemap \
    && conda install -c conda-forge cartopy

ENV PROJ_LIB=/usr/share/proj

RUN /opt/conda/bin/pip --disable-pip-version-check --no-cache-dir install requests matplotlib numpy scipy aacgmv2 image

# make necessary directories
RUN mkdir Output
RUN chown ovation:ovation Output

USER $USERNAME

# Add the source code into the image
COPY --chown=ovation:ovation source source
COPY --chown=ovation:ovation source_plots source_plots
COPY --chown=ovation:ovation driver.sh .

COPY --chown=ovation:ovation SW_Data SW_Data
COPY --chown=ovation:ovation configuration configuration
