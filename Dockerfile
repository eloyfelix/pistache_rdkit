FROM debian:bookworm-slim

ARG RDKIT_VERSION=Release_2024_03_4
ARG PISTACHE_COMMIT=46ffbb8228333b146c0ac2294439ba68623890b3

RUN apt-get update --fix-missing && \
    apt-get install -y g++ \
                       cmake \
                       pkg-config \
                       meson \
                       python3-setuptools \
                       libssl-dev \
                       curl \
                       unzip \
                       git \
                       zlib1g-dev \
                       nlohmann-json3-dev \
                       libboost-dev \
                       libboost-iostreams-dev \
                       libboost-system-dev \
                       libboost-serialization-dev && \
    apt-get -qq -y autoremove && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

# Install RDKit
RUN curl -LO https://github.com/rdkit/rdkit/archive/${RDKIT_VERSION}.tar.gz && \
    tar -xzf ${RDKIT_VERSION}.tar.gz && \
    mv rdkit-${RDKIT_VERSION} rdkit && \
    mkdir rdkit/build && \
    cd /rdkit/build && \
    cmake -Wno-dev \
          -DCMAKE_BUILD_TYPE=Release \
          -DRDK_BUILD_INCHI_SUPPORT=ON \
          -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
          -DCMAKE_INSTALL_PREFIX=/usr \
          -DCMAKE_SYSTEM_PREFIX_PATH=/usr \
          -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
          -DRDK_INSTALL_INTREE=OFF \
          -DRDK_BUILD_CPP_TESTS=OFF \
          .. && \
    make -j $(nproc) && \
    make install && \
    rm -rf /rdkit /${RDKIT_VERSION}.tar.gz

# Install pistache
RUN git clone https://github.com/pistacheio/pistache.git && \
    cd pistache && \
    git checkout ${PISTACHE_COMMIT} && \
    meson setup build \
      --buildtype=release \
      -DPISTACHE_USE_SSL=true \
      -DPISTACHE_BUILD_EXAMPLES=false \
      -DPISTACHE_BUILD_TESTS=false \
      -DPISTACHE_BUILD_DOCS=false \
      --prefix=/usr && \
    meson compile -C build && \
    meson install -C build

COPY src app/src
COPY CMakeLists.txt app/CMakeLists.txt

RUN mkdir app/build && \
    cd app/build && \
    cmake .. && \
    make && \
    cp PistacheRDKit ../../PistacheRDKit

CMD ["./PistacheRDKit", "9080", "4"]
