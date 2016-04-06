# pipelign

Pipeline for aligning virus sequences

Typing ``python3 Pipelign.py -h`` will give a short usage information

Dependencies:

- [Python 3](https://www.python.org/)
- [BioPython](http://biopython.org/wiki/Main_Page): for sequence processing
- [CD-HIT](http://weizhongli-lab.org/cd-hit/): for clustering 'long' sequences
- [IQ-TREE](http://www.cibiv.at/software/iqtree/): for generating a phylogenetic tree of representative 'long' sequences
- [MAFFT](http://mafft.cbrc.jp/alignment/software/): for alignment
- [HMMER](http://hmmer.org/): for matching sequence fragments to clusters
