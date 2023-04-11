# Minimal Docker image for MAPLE v0.3.1 using Ubuntu base
FROM ubuntu:20.04

# install MAPLE
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y unzip wget zip && \

    # install pypy
    wget -qO- "https://downloads.python.org/pypy/pypy3.9-v7.3.11-linux64.tar.bz2" | tar -xj && \
    mv pypy*/bin/* /usr/local/bin/ && \
    mv pypy*/include/* /usr/local/include/ && \
    mv pypy*/lib/* /usr/local/lib/ && \
    rm -rf pypy* && \

    # install MAPLE
    wget -q "https://github.com/NicolaDM/MAPLE/archive/refs/heads/main.zip" && \
    unzip main.zip && \
    rm main.zip && \
    mv MAPLE-main /usr/local/bin/MAPLE && \
    sed -i '1s/^/#! \/usr\/bin\/env pypy3\n/' /usr/local/bin/MAPLE/MAPLEv*.py && \
    chmod a+x /usr/local/bin/MAPLE/MAPLEv*.py && \
    ln -s /usr/local/bin/MAPLE/MAPLEv*.py /usr/local/bin/MAPLE.py && \

    # clean up
    rm -rf /root/.cache /tmp/*
