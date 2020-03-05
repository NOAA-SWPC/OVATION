FROM python:3.8-slim-buster

ARG USERNAME=ovation
ARG USER_UID=1001
ARG USER_GID=$USER_UID

# Add ovation user/group
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME

RUN pip --disable-pip-version-check --no-cache-dir install requests
RUN pip --disable-pip-version-check --no-cache-dir install matplotlib numpy scipy
USER $USERNAME

# Add the source doe into the image
COPY --chown=ovation:ovation . /
