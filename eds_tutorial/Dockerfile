FROM whitelab/plumed-gromacs
USER root
RUN apt-get update; exit 0;
RUN apt-get install -y python3-numpy python3-scipy python3-pip python3-tk
RUN pip3 install matplotlib numpy 
USER whitelab
WORKDIR $HOME/scratch
