FROM ubuntu:20.04

MAINTAINER Yunpeng Men "ypmen@mpifr-bonn.mpg.de"

SHELL ["/bin/bash", "-c"] 

RUN useradd -ms /bin/bash pulsarx

ENV DEBIAN_FRONTEND noninteractive
ENV HOME /home/pulsarx

ENV LD_LIBRARY_PATH=$HOME/software/lib:$LD_LIBRARY_PATH
ENV YMW16_DIR=$HOME/software/PulsarX/src/ymw16
ENV PATH=$PATH:$HOME/software/bin
ENV OMP_NUM_THREADS=1
ENV PSRCAT_FILE=/home/pulsarx/software/psrcat_tar/psrcat.db
ENV TEMPO2=/home/pulsarx/software/tempo2/T2runtime

USER root

RUN apt-get update

RUN apt-get install -y git \
    autoconf \
    libtool \
    g++ \
    gfortran \
    libboost-program-options-dev \
    libboost-log-dev \
    libcfitsio-dev \
    libblas-dev \
    liblapack-dev \
    libeigen3-dev \
    libpng-dev \
    libfftw3-dev \
    pgplot5 \
    wget \
    make \
    libx11-dev \
    python3-dev \
    python3-distutils \
    python3-numpy \
    python3-matplotlib

USER pulsarx

#install sofa
WORKDIR $HOME
RUN mkdir software
WORKDIR $HOME/software
RUN wget https://www.iausofa.org/2020_0721_C/sofa_c-20200721.tar.gz --no-check-certificate && \
    tar -zxvf sofa_c-20200721.tar.gz
WORKDIR $HOME/software/sofa/20200721/c/src
RUN make && make test

#install tempo2
WORKDIR $HOME/software
RUN git clone https://bitbucket.org/psrsoft/tempo2.git
WORKDIR $HOME/software/tempo2
RUN rm -rf .git
RUN ./bootstrap
RUN ./configure --prefix=$HOME/software
RUN make && make install
RUN make plugins-install
RUN make clean

#install psrcat
WORKDIR $HOME/software
RUN wget https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz && \
    tar -zxvf psrcat_pkg.tar.gz
WORKDIR $HOME/software/psrcat_tar
RUN source makeit && cp psrcat $HOME/software/bin

WORKDIR $HOME/software
#install PlotX
RUN git clone https://github.com/ypmen/PlotX.git
#install PulsarX
RUN git clone https://github.com/ypmen/PulsarX.git
#install BasebandX
RUN git clone https://github.com/ypmen/BasebandX.git
#install TransientX
RUN git clone https://github.com/ypmen/TransientX

WORKDIR $HOME/software/PlotX
RUN ./bootstrap
RUN ./configure --prefix=$HOME/software
RUN make && make install
RUN make clean

WORKDIR $HOME/software/PulsarX
RUN ./bootstrap
RUN ./configure --prefix=$HOME/software CXXFLAGS="-std=c++11 -O3" LDFLAGS="-L$HOME/software/sofa/20200721/c/src -L$HOME/software/lib" CPPFLAGS="-I$HOME/software/sofa/20200721/c/src -I$HOME/software/include"
RUN make && make install
RUN make clean

WORKDIR $HOME/software/BasebandX
RUN ./bootstrap
RUN ./configure --prefix=$HOME/software CXXFLAGS="-std=c++11 -O3"
RUN make && make install
RUN make clean

WORKDIR $HOME/software/TransientX
RUN ./bootstrap
RUN ./configure --prefix=$HOME/software CXXFLAGS="-std=c++11 -O3" LDFLAGS="-L$HOME/software/sofa/20200721/c/src -L$HOME/software/lib" CPPFLAGS="-I$HOME/software/sofa/20200721/c/src -I$HOME/software/include"
RUN make && make install
RUN make clean

WORKDIR $HOME
