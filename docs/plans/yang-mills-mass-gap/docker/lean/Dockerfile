# Lean image for the Yang-Mills audit. Pins the toolchain of the proof under review.
# Adapted from claw-work/Scripts/lean4-docker (fixes ELAN_HOME, which Docker does not expand from $HOME).
# Build:  PLAN/bin/ym-lean build
FROM debian:bookworm-slim

ARG LEAN_TOOLCHAIN=leanprover/lean4:v4.32.0
ARG UID=1000
ARG GID=1000

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl git build-essential xz-utils zstd python3 \
 && rm -rf /var/lib/apt/lists/*

RUN groupadd -g ${GID} lean && useradd -m -u ${UID} -g ${GID} -s /bin/bash lean
USER lean
ENV ELAN_HOME=/home/lean/.elan
ENV PATH=/home/lean/.elan/bin:${PATH}

RUN curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o /tmp/elan-init.sh \
 && sh /tmp/elan-init.sh -y --default-toolchain none \
 && rm /tmp/elan-init.sh \
 && elan toolchain install ${LEAN_TOOLCHAIN} \
 && elan default ${LEAN_TOOLCHAIN}

WORKDIR /work
CMD ["lean", "--version"]
