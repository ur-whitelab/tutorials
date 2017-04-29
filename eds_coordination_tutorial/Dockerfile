FROM whitelab/plumed-lammps
USER root

RUN apt-get update; exit 0;
RUN apt-get install -y python3-numpy python3-scipy python3-pip python3-tk
RUN pip3 install matplotlib

USER whitelab
WORKDIR $HOME/scratch
RUN wget https://gist.githubusercontent.com/whitead/ba45b7aff24bf353ae93a96dfc3855d2/raw/ae39891e85a5573eae5e7605576616410495867f/convert_gr.py
