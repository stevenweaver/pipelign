import luigi
import sys, os, shutil, subprocess, argparse, textwrap
import utils
from Bio import SeqIO, Phylo


class runCDHIT(luigi.Task):

  '''
    CD-HIT is used to group similar sequences together in clusters for alignment
  '''
  task_namespace = 'main'

  longName = luigi.Parameter()
  alphabet = luigi.Parameter()
  per = luigi.IntParameter()
  thread = luigi.IntParameter()
  cDir = luigi.Parameter()
  tName = luigi.Parameter()
  zName = luigi.Parameter()

  def output(self):
    return luigi.LocalTarget("grp")

  def run(self):

    print('{task} says: Running CD-HIT/CD-HIT-EST to group long sequences into clusters based on sequence similarity')

    seqs = list(SeqIO.parse(self.longName,'fasta'))

    if len(seqs) < 2:
      try:
        shutil.copy(self.longName,'grp')
      except OSError as e:
        print(e)
        utils.cZip(self.cDir,self.tName,self.zName)
      return


    # create log file for cdhit
    lh = open('cdhit.log','w')

    if self.alphabet == 'dna' or self.alphabet == 'rna':
      cl = 'cd-hit-est -c %f -n 5 -i %s -o grp -d 0 -T %d' % (self.per,self.longName,self.thread)

    elif alphabet == 'aa':
      cl = 'cd-hit -c %f -n 5 -i %s -o grp -d 0 -T %d' % (self.per,self.longName,self.thread)

    print("{task} says: executing command " + cl)

    try:
      subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
      print(e)
      utils.cZip(self.cDir,self.tName,self.zName)

    lh.close() # close log file

    print('{task} says: CD-HIT created two files:')
    print('{task} says: <grp> for cluster representative long sequences')
    print('{task} says: <grp.clstr> for cluster assignment of long sequences')


#*************************************************************************
class makeClusters(luigi.Task):

  '''
    Separate files are created for each clusters
  '''

  task_namespace = 'main'

  longName = luigi.Parameter()

  def output(self):
    return luigi.LocalTarget('clusterList.txt')

  def requires():
    return runCDHIT(longName=longName, alphabet=alphabet, per=per, thread=thread, cDir=cDir, tName=tName, zName=zName)

  def run(self):

    print('\nCreating cluster files for long sequences')
    seqs = list(SeqIO.parse(self.longName,'fasta')) # load all the full sequences
    clsSize = list()

    if len(seqs) < 2:
      shutil.copy(self.longName,'grp.0.fas')
      clsSize.append(1)
      fh = open('clusterList.txt','w')
      fh.write('%s\t0' %seqs[0].id)
      fh.close()
      return 1, clsSize

    cName = 'grp.clstr'
    # read cluster file
    lines = [line.strip('\n') for line in open(cName,'r')]
    # flag for the beginning of first cluster list
    start = 0
    # hold sequences of a cluster
    cSeq = []
    # count clusters
    cls = 0
    # clusterList string
    st = ''
    # IDs for the full sequences
    ids = []

    for seq in seqs:
      ids.append(seq.id)

    '''
      read cluster file and make separate cluster files
    '''

    if 'Cluster' not in lines[0]:
      msg = '\n\n"grp.clstr" does not contain any cluster'
      msg += '\nPlease try running the program again\n'
      sys.exit(msg)

    for i in range(1, len(lines)):
      if 'Cluster' in lines[i]: # start of a new cluster list
        gName = 'grp.' + str(cls) + '.fas'
        cls = cls + 1
        SeqIO.write(cSeq,gName,'fasta')
        cSeq = []
        clsSize.append(len(cSeq))
      else: # continue with the existing cluster
        seqID = lines[i].split()[2].replace('>','').replace('...','')
        st += seqID + '\t' + str(cls) + '\n' # updating clusterList file content
        sInd = ids.index(seqID)
        cSeq.append(seqs[sInd])

    gName = 'grp.' + str(cls) + '.fas'
    cls = cls + 1
    SeqIO.write(cSeq,gName,'fasta')

    fh = open('clusterList.txt','w')
    fh.write(st)
    fh.close()

    print('\tTotal %d cluster file(s) created; example name <cls.0.fas>' % cls)

    return cls, clsSize

#*************************************************************************
class addClusterNumberToReps(luigi.Task):

  '''
    - Reads in the cluster representative FASTA file and the clusterList.txt file
    - Adds cluster number to the sequence header e.g. >seq1_0
    - Temporary file <clsReps.fas> is written
  '''

  task_namespace = 'main'

  repName = luigi.Parameter()
  lstFile = luigi.Parameter()
  outFile = luigi.Parameter()
  clsSize = luigi.IntParameter()

  def requires(self):
    return makeClusters(longName=longName)

  def output(self):
    return luigi.LocalTarget(self.outFile)

  def run(self):

    print('\nWriting cluster representatives with cluster number in the header')

    cList = [line.strip() for line in open(self.lstFile,'r')]

    cID = [] # sequence ids
    cNum = [] # cluster numbers

    for line in cList:
      words = line.split()
      cID.append(words[0])
      cNum.append(words[1])

    seqs = list(SeqIO.parse(self.repName,'fasta'))

    for seq in seqs:
      if seq.id in cID:
        ind = cID.index(seq.id)
        seq.id = seq.id + '_' + str(cNum.count(cNum[ind])) + '_' + cNum[ind]
        seq.name = seq.id
        seq.description = seq.id
      else:
        sys.exit('\nSequence %s does not have a cluster. Pipelign is exiting' % seq.id)

    SeqIO.write(seqs,self.outFile,'fasta')



if __name__ == '__main__':
  pass
  # luigi.run(['main.runCDHIT', '--local-scheduler', '--longName ./tests/rsrc/input.fas'])

