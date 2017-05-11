FROM whitelab/plumed-lammps
USER root

RUN apt-get update; exit 0;
RUN apt-get install -y python3-numpy python3-scipy python3-pip python3-tk wget && apt-get clean;
RUN pip3 install --no-cache-dir matplotlib numpy MDAnalysis



USER whitelab
WORKDIR $HOME/scratch

