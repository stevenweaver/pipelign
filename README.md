# Pipelign

Pipeline for aligning virus sequences

Typing ``python3 Pipelign.py -h`` will give a short usage information

Dependencies:

- [Python 3](https://www.python.org/)
- [BioPython](http://biopython.org/wiki/Main_Page): for sequence processing
- [CD-HIT](http://weizhongli-lab.org/cd-hit/): for clustering 'long' sequences
- [IQ-TREE](http://www.cibiv.at/software/iqtree/): for generating a phylogenetic tree of representative 'long' sequences
- [MAFFT](http://mafft.cbrc.jp/alignment/software/): for alignment
- [HMMER](http://hmmer.org/): for matching sequence fragments to clusters

## Docker

Docker repository for **Pipelign** is hosted on [asmmhossain/Pipelign](https://hub.docker.com/r/asmmhossain/pipelign/). Anyone having [Docker](https://www.docker.com/) installed can get **Pipelign** using the command:

``docker pull asmmhossain/pipelign``

The command for launching **Pipelign** in **docker** is:

``docker run -i -t --rm -v $(pwd):/data asmmhossain/pipelign /sbin/my_init -- bash -l``

Please make sure the sequence file resides in the current directory from where the command is executed. Output file will be stored in the current directory.


