FROM ubuntu

RUN apt-get update
RUN apt-get install -y \
        git \
        gcc \
        make \
        libgfortran-4.8-dev \
        libblas-dev \
        liblapack-dev \
        python2.7

RUN ln -s $(which python2.7) $(dirname $(which python2.7))/python
RUN git clone https://github.com/dandavison/shellfish.git
WORKDIR shellfish/src
RUN make all
WORKDIR ..
