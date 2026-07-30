FROM python:3.10-slim-bookworm AS builder

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        cmake \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cmake -S rmath -B build-rmath \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/rmath \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    && cmake --build build-rmath --parallel \
    && cmake --install build-rmath \
    && cmake -S . -B build-container \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH=/opt/rmath \
        -DPANGROWTH_WITH_GGCAT=OFF \
    && cmake --build build-container --parallel


FROM python:3.10-slim-bookworm AS runner

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        libgomp1 \
        tar \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /opt/pangrowth/requirements.txt
RUN pip install --no-cache-dir --requirement /opt/pangrowth/requirements.txt

COPY --from=builder /src/build-container/pangrowth /usr/local/bin/pangrowth
COPY scripts/plot_core.py scripts/plot_growth.py scripts/plot_hist.py /opt/pangrowth/scripts/
RUN chmod 0755 /usr/local/bin/pangrowth /opt/pangrowth/scripts/*.py

ENV PATH="/opt/pangrowth/scripts:${PATH}" \
    MPLCONFIGDIR="/tmp/matplotlib"

WORKDIR /work

# CloWM/Nextflow supplies the command, so this image intentionally has no
# ENTRYPOINT. The harmless Python base-image CMD is overridden by Nextflow.
