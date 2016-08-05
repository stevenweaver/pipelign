import luigi
import sys, os, shutil, subprocess, argparse, textwrap
import utils
from Bio import SeqIO, Phylo

class makeIQTree():
  '''
    - Constructs phylogenetic tree using IQ-TREE
  '''

  task_namespace = 'main'

  alnFile = luigi.Parameter()
  thread = luigi.IntParameter()
  cDir = luigi.Parameter()
  tName = luigi.Parameter()
  zName = luigi.Parameter()

  def requires():
    return makeClusters(longName=longName)

  def output():
    return addClusterNumberToReps(longName=longName)

  def run():
    lh = open('iqtree.log','a')

    print('\nCreating IQ-TREE from %s' % alnFile)

    cl = 'iqtree-omp -s %s -m TEST -nt %s' % (alnFile,thread)

    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      print(e)
      cZip(cDir,tName,zName)

    lh.close()

    print('\tTree file is written in %s.treefile' % alnFile)
    print('\tLog file is written in %s.log' % alnFile)

#***********************************************************************

#***********************************************************************
def drawMidPointRootTree(treeFile):
  '''
   - Displays a dendogram of the tree generated from cluster representatives
  '''
  task_namespace = 'main'

  alnFile = luigi.Parameter()
  thread = luigi.IntParameter()
  cDir = luigi.Parameter()
  tName = luigi.Parameter()
  zName = luigi.Parameter()


  def run():
    print('\nThe phylogenetic tree for the cluster representatives is shown below:\n')
    tree = Phylo.read(treeFile,'newick')
    tree.root_at_midpoint()
    Phylo.draw_ascii(tree)
    print('\n')
    '''
    try:
      input('press ENTER to continue: ')
    except SyntaxError:
      pass
    '''


if __name__ == '__main__':
  pass
  # luigi.run(['main.runCDHIT', '--local-scheduler', '--longName ./tests/rsrc/input.fas'])

