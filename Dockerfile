FROM python:3.8-slim-buster

ARG USERNAME=ovation
ARG USER_UID=1001
ARG USER_GID=$USER_UID

# Add ovation user/group
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME

RUN pip --disable-pip-version-check --no-cache-dir install requests
RUN pip --disable-pip-version-check --no-cache-dir install matplotlib numpy scipy

# make necessary directories
#RUN mkdir SW_Data && \
#  mkdir configuration && \
#  mkdir source_plots && \
RUN mkdir Output
RUN chown ovation:ovation Output
#  mkdir source

USER $USERNAME

# Add the source code into the image
COPY --chown=ovation:ovation source source

COPY --chown=ovation:ovation SW_Data SW_Data
COPY --chown=ovation:ovation configuration configuration
COPY --chown=ovation:ovation source_plots source_plots

WORKDIR source
