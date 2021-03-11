FROM ubuntu:18.04
MAINTAINER Jaswinder Singh (jaswinder.singh3@griffithuni.edu.au)

RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN apt-get update && apt-get install -y build-essential wget virtualenv git python-minimal cpanminus gawk
RUN cpanm Graph

RUN wget 'https://www.dropbox.com/s/p94grd6c0v1eg73/SPOT-RNA2.tar.xz' || wget 'https://app.nihaocloud.com/f/3e826caf8efc43adaaa0/?dl=1' && tar -xvf SPOT-RNA2.tar.xz && rm SPOT-RNA2.tar.xz
WORKDIR SPOT-RNA2

RUN wget -O utils/models_ckps.tar.xz 'https://www.dropbox.com/s/udzcsva76lh5wvq/models_ckps.tar.xz' || wget -O utils/models_ckps.tar.xz 'https://app.nihaocloud.com/f/586acb2658d74ccb92b8/?dl=1' && tar -xvf utils/models_ckps.tar.xz -C utils/ && rm utils/models_ckps.tar.xz
RUN virtualenv -p python3.6 venv && source ./venv/bin/activate &&  pip install tensorflow==1.14.0 && pip install -r requirements.txt && deactivate

RUN wget 'eddylab.org/infernal/infernal-1.1.3-linux-intel-gcc.tar.gz' && tar -xvzf infernal-*.tar.gz && rm infernal-*.tar.gz
RUN wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*+-x64-linux.tar.gz' && tar -xvzf ncbi-blast-*+-x64-linux.tar.gz && rm ncbi-blast-*+-x64-linux.tar.gz
RUN git clone https://github.com/jaswindersingh2/SPOT-RNA.git && cd SPOT-RNA && wget 'https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz' || wget -O SPOT-RNA-models.tar.gz 'https://app.nihaocloud.com/f/fbf3315a91d542c0bdc2/?dl=1' && tar -xvzf SPOT-RNA-models.tar.gz && rm SPOT-RNA-models.tar.gz && cd ../
RUN git clone "https://github.com/sokrypton/GREMLIN_CPP" && cd GREMLIN_CPP && g++ -O3 -std=c++0x -o gremlin_cpp gremlin_cpp.cpp -fopenmp && cd ../
RUN git clone 'https://github.com/LinearFold/LinearPartition.git' && cd LinearPartition/ && make && cd ../
