#!/usr/bin/env python3
'''
  Pipelign.py

  A python based program to align virus sequences.

  The program takes as input a single FASTA formatted file

  Returns a FASTA formatted alignment file

'''
#*********************************************************************

import sys, os, shutil, subprocess, argparse, textwrap
from Bio import SeqIO
import tempfile
import time
from Bio import Phylo

#************************************
class Struct:
  '''
   - this class creates a structure of the parameter values needed for different functions.
   - parameters can be accessed as structName.itemName e.g. mArgs.gCode 
  '''
  def __init__(self, **kwargs):
    for k, v in kwargs.items():
      setattr(self,k,v)
      
#*************************************      

#*************************************
class MyStruct(Struct):
  pass

#*****************************************

#******************************************
def lengthThreshold(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % (x,))
  return x
#***************************************

#******************************************
def iterThreshold(x):
  x = int(x)
  if x < 1:
    raise argparse.ArgumentTypeError('%r must be a non-zero positive integer' % (x,))
  return x
#***************************************


#****************************************
def genCodeLimit(x):
  x = int(x)
  gCodes = list(1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25)
  if x not in gCodes:
    raise argparse.ArgumentTypeError('%r is not a valid genetic code' % (x,))
  else:
    return x

#****************************************

def cZip(cDir,tName,zName):
  '''
  creates a zip file of the temporary directory
  '''  
  os.chdir(cDir)
  
  #zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S') 
  try:
    shutil.make_archive(zName,'zip',tName)
  except: 
    sys.exit(e)
  
  print('\nArchive for all temporary files created in %s.zip\n' % zName)
  sys.exit()

#******************************************

#***********************************************************************

def ginsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName):
  '''
    - aligns sequences using MAFFT's G-INS-i method
  '''
  
  lh = open(log,'a')
  
  cl = 'mafft --globalpair --thread %d --maxiterate %d --preservecase %s > %s' % (thread, mIterLong,seqFile, alnFile)

  try:
    subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)

  lh.close()
#***********************************************************************

#***********************************************************************

def linsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName):
  '''
    - aligns sequences using MAFFT's L-INS-i method
  '''
  
  lh = open(log,'a')
  
  cl = 'mafft --localpair --thread %d --maxiterate %d --preservecase %s > %s' % (thread, mIterLong,seqFile, alnFile)

  try:
    subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
    
  lh.close()
#***********************************************************************

#***********************************************************************

def fftnsi(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName):
  '''
    - aligns sequences using MAFFT's FFT-NS-i method
  '''
  
  lh = open(log,'a')
  
  cl = 'mafft --thread %d --maxiterate %d --preservecase %s > %s' % (thread, mIterLong,seqFile, alnFile)

  try:
    subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)

  lh.close()
#***********************************************************************

#***********************************************************************

def fftns2(seqFile,alnFile,thread,mIterLong,log,cDir,tName,zName):
  '''
    - aligns sequences using MAFFT's FFT-NS-i method
  '''
  
  lh = open(log,'a')
  
  cl = 'mafft --thread %d --retree 2 --preservecase %s > %s' % (thread, seqFile, alnFile)

  try:
    subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  lh.close()
#***********************************************************************
#***********************************************************************

def addFragments(fName,aName,oName,thread,log,cDir,tName,zName):
  '''
    - add fragments to the alignment using MAFFT's --add
  '''

  lh = open(log,'a')
  
  cl = 'mafft --preservecase --thread %d --maxiterate 1000 --addfragments %s %s > %s' % (thread, fName, aName, oName)
        
  try:
    subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  lh.close()
#***********************************************************************

#************************************************************************

def deAlign(iFile, dFile):
  '''
  Removes gaps (if any) from the input sequence file
  
  ''' 
  
  print("\nRemoving the gap characters '-'/'.' (if any) from the sequences")
  seqs = list(SeqIO.parse(iFile,'fasta'))
  
  st = ''
  
  for seq in seqs:
    st += '>' + seq.id + '\n' + str(seq.seq).replace('-','').replace('.','') + '\n'
  
  fh = open(dFile,'w')
  fh.write(st)
  fh.close()
  
  print('\tGapless sequence file written in <%s>' % dFile)
  
#*************************************************************************

#*************************************************************************
def getLengths(inFile):
  '''
  '''

#************************************************************************

def separateFullFragment(iFile, thr, longName, fragName):
  '''
    Reads in the input sequence file.
      - finds length of the longest sequence
      - find minimum length required for full length sequences
      - writes full length sequences and fragments into two separate files 
  '''
  fEmpty = 1
  
  print('\nCreating separate files for long sequences and fragments')
  seqs = list(SeqIO.parse(iFile,'fasta'))
  
  mLen = -1
  
  for seq in seqs:
    if len(seq.seq) > mLen:
      mLen = len(seq.seq)
  
  minLengthFull = int(thr * mLen)
  
  full = []
  frag = []
  
  for seq in seqs:
    if len(seq.seq) < minLengthFull:
      frag.append(seq)
    else:
      full.append(seq)
      
  
  SeqIO.write(full,longName,'fasta') 
  print('\tLong sequences are written in <%s>' % longName)

  if len(frag) > 0:
    SeqIO.write(frag, fragName,'fasta')
    print('\tFragments are written in <%s>' % fragName)
    fEmpty = 0
    
  
  if fEmpty:
    return 1
  else:
    return 0
#************************************************************************

#************************************************************************
def runCDHIT(longName,alphabet,per,thread,cDir,tName,zName):
  '''
    CD-HIT is used to group similar sequences together in clusters for alignment
    
  '''
   
  seqs = list(SeqIO.parse(longName,'fasta'))
  
  if len(seqs) < 2:
    try:
      shutil.copy(longName,'grp')
    except OSError as e:
      print(e)
      cZip(cDir,tName,zName)
    return
  
  print('\nRunning CD-HIT/CD-HIT-EST to group long sequences into clusters based on sequence similarity')
  
  lh = open('cdhit.log','w') # create log file for cdhit
  
  if alphabet == 'dna' or alphabet == 'rna':
    cl = 'cd-hit-est -c %f -n 5 -i %s -o grp -d 0 -T %d' % (per,longName,thread)
    #print(cl)

  elif alphabet == 'aa':
    cl = 'cd-hit -c %f -n 5 -i %s -o grp -d 0 -T %d' % (per,longName,thread)
    #print(cl)
    
  try:
    subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  lh.close() # close log file   
  
  print('\tCD-HIT created two files:')
  print('\t\t  <grp> for cluster representative long sequences')
  print('\t\t  <grp.clstr> for cluster assignment of long sequences')
#*************************************************************************

#*************************************************************************
def makeClusters(longName):
  '''
    Separate files are created for each clusters
  '''
  
  print('\nCreating cluster files for long sequences')
  
  seqs = list(SeqIO.parse(longName,'fasta')) # load all the full sequences

  clsSize = list()
  
  if len(seqs) < 2:
    shutil.copy(longName,'grp.0.fas')
    clsSize.append(1)
    fh = open('clusterList.txt','w')
    fh.write('%s\t0' %seqs[0].id)
    fh.close()
    return 1, clsSize
  
  cName = 'grp.clstr'
  
  lines = [line.strip('\n') for line in open(cName,'r')] # read cluster file
  
  start = 0 # flag for the beginning of first cluster list
   
  cSeq = [] # hold sequences of a cluster
  
  cls = 0 # count clusters
  
  st = '' # clusterList string
  
  ids = [] # IDs for the full sequences
  
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
#***********************************************************************

#***********************************************************************

def addClusterNumberToReps(repName,lstFile,outFile,clsSize):
  '''
    - Reads in the cluster representative FASTA file and the clusterList.txt file
    - Adds cluster number to the sequence header e.g. >seq1_0  
    - Temporary file <clsReps.fas> is written
  '''
  
  print('\nWriting cluster representatives with cluster number in the header')
  
  cList = [line.strip() for line in open(lstFile,'r')]
  
  cID = [] # sequence ids 
  cNum = [] # cluster numbers
  
  for line in cList:
    words = line.split()
    cID.append(words[0])
    cNum.append(words[1])
    
  seqs = list(SeqIO.parse(repName,'fasta'))
  
  for seq in seqs:
    if seq.id in cID:
      ind = cID.index(seq.id)
      
      seq.id = seq.id + '_' + str(cNum.count(cNum[ind])) + '_' + cNum[ind] 
      seq.name = seq.id
      seq.description = seq.id
    else:
      sys.exit('\nSequence %s does not have a cluster. Pipelign is exiting' % seq.id)
  
  SeqIO.write(seqs,outFile,'fasta')    

#***********************************************************************

#***********************************************************************

def makeClusterRepsAlignment(repFile,outFile,thread,mIterL,cDir,tName,zName):
  '''
    - Creates a multiple sequence alignment from cluster reps with cluster numbers
    - uses G-INS-i | L-INS-i | FFT-NS-i 
    - writes the alignment in <outFile> 
  '''

  print('\nAligning cluster representative sequences')
  
  #ginsi(repFile,outFile,thread,mIterL,'clsRepAln.log',cDir,tName,zName)
  fftnsi(repFile,outFile,thread,mIterL,'clsRepAln.log',cDir,tName,zName)
  
  print('\tAlignment written in %s' % outFile)
#***********************************************************************

#***********************************************************************

def makeIQTree(alnFile,thread,cDir,tName,zName):
  '''
    - Constructs phylogenetic tree using IQ-TREE
  '''
  
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
def drawTree(treeFile):
  '''
   - Displays a dendogram of the tree generated from cluster representatives
  '''
  
  print('\nThe phylogenetic tree for the cluster representatives is shown below:\n')
  tree = Phylo.read(treeFile,'newick')
  Phylo.draw_ascii(tree)
  print('\n')
  '''
  try:
    input('press ENTER to continue: ')
  except SyntaxError:
    pass
  '''
#***********************************************************************

#***********************************************************************
def drawMidPointRootTree(treeFile):
  '''
   - Displays a dendogram of the tree generated from cluster representatives
  '''
  
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
#***********************************************************************

#***********************************************************************
def alnFullSequenceClusters(nClusters,thread,mIterL,cDir,tName,zName):
  '''
    Full sequences in each clusters will be aligned using L-INS-i/clustalo
  
  '''
  log = 'clusterAlign.log'
  lh = open(log,'w') # open log file for cluster alignment
  lh.close()
  
  print('\nAligning clusters')
  for i in range(nClusters):
    cName = 'grp.' + str(i) + '.fas'
    aName = 'grp.' + str(i) + '.aln'
    seqs = list(SeqIO.parse(cName,'fasta'))
    
    print('\tAligning cluster %d of %d sequences' % (i,len(seqs)))
    
    if len(seqs) > 1:
      #print('\tAligning cluster %d of %d sequences' % (i,len(seqs)))
      #ginsi(cName,aName,thread,mIterL,log,cDir,tName,zName)
      #linsi(cName,aName,thread,mIterL,log,cDir,tName,zName)
      fftnsi(cName,aName,thread,mIterL,log,cDir,tName,zName)
    else:
      shutil.copyfile(cName,aName)

  #lh.close() # close log file
#***********************************************************************  

#***********************************************************************  
def makeHMMdb(nClusters,cDir,tName,zName,thread,lFile,alpha):
  '''
    Create profile HMMs for each of the full length alignments 
    Create HMM_DB from the pHMMs
     
  '''
  #if fragEmpty:
    #return # if there is no fragment sequences

  
  print('\nCreating profile HMMs for %d cluster alignments' % nClusters)
  
  lh = open(lFile,'a') # create log file for hmmer
  
  for i in range(nClusters):
    aName = 'grp.' + str(i) + '.aln'
    hName = 'grp.' + str(i) + '.hmm'
    
    #*****testing for failure of hmmbuild
    # hmmbuild fails when alignment does not have any gaps
    # solution: add a '-' at the end of each sequences ????
    aseqs = list(SeqIO.parse(aName,'fasta'))

    if len(aseqs) == 1:
      aseqs.append(aseqs[0])
      aseqs[1].id = aseqs[0].id + '_1'
      aseqs[1].name = aseqs[1].id
      aseqs[1].description = aseqs[1].id

    if len(aseqs) > 1:
      sumGap = 0
      for seq in aseqs:
        sumGap = sumGap + seq.seq.count('-')
      if sumGap == 0:
        for seq in aseqs:
          seq.seq = seq.seq + '-'
        aName = 'grp.' + str(i) + '.tmp.aln'
        SeqIO.write(aseqs,aName,'fasta')  
        #aName = 'temp.aln'
    
    #**************
    if alpha == 'dna':
      cl = 'hmmbuild --dna --cpu %d %s %s' % (thread, hName, aName)
    
    if alpha == 'rna':
      cl = 'hmmbuild --rna --cpu %d %s %s' % (thread, hName, aName)
    
    elif alpha == 'aa':
      cl = 'hmmbuild --amino --cpu %d %s %s' % (thread, hName, aName)
    
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      #sys.exit(e)
      print(e)
      cZip(cDir,tName,zName)
    print('\t<%s> created' % hName)
  
  print('\nCreating the HMM database from the profile HMMs')
  
  dName = 'grp.hmm'
  
  fh = open(dName,'w')
  for i in range(nClusters):
    hName = 'grp.' + str(i) + '.hmm'
    shutil.copyfileobj(open(hName,'r'),fh)
  fh.close()
  
  # create HMM database
  
  cl = 'hmmpress %s' % dName
  
  try:
    subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    #sys.exit(e)
    print(e)
    cZip(cDir,tName)

  print('\tHMM database <%s> is created' % dName)
    
  lh.close() # close log file     
#***********************************************************************

#***********************************************************************  
def searchHMMdb(lFile,thread,fragEmpty,alpha,res,cDir,tName,zName):
  '''
    HMM database is searched with the fragments to assign them a cliuster for alignment
  '''  
  lh = open(lFile,'a')
  
  if not fragEmpty:
    print('\nFragments will be searched against the HMM database to assign clusters')
    print('\tResults of the database search are being written on <%s>' % res)
    
    if alpha == 'dna' or alpha == 'rna':
      cl = 'nhmmscan --cpu %d --tblout %s --noali grp.hmm frag.fas' % (thread,res)
      
    elif alpha == 'aa':
      cl = 'hmmscan --cpu %d --tblout %s --noali grp.hmm frag.fas' % (thread,res)
    
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      print(e)
      cZip(cDir,tName,zName)
  
  lh.close()
#**********************************************************************

#**********************************************************************    
def parseHMMsearchResult(nClusters,fragFile,res,keepOrphans):
  '''
    Reads in the result file from hmm.out file to determine cluster for fragments
  '''
  
  #global addNClusters
  print('\nWriting cluster files for fragments')
  fSeqIDs = [] # IDs of the fragments
  
  fFlag = 0
  
  #if fragEmpty:
    #return # if there is no fragment sequences
        
  fSeqs = list(SeqIO.parse(fragFile,'fasta'))
  
  for seq in fSeqs:
    fSeqIDs.append(seq.id)
  
  fIndex = [-1 for x in range(len(fSeqIDs))] # to hold cluster number assigned for each fragments
  
  lines = [line.strip('\n') for line in open(res,'r')]
  
  for l in lines:
    if 'grp' in l and l[0] != '#':
      words = l.split()
      g = words[0].split('.')[1]
      q = words[2]

      ind = fSeqIDs.index(q)
      if fIndex[ind] == -1:
        fIndex[ind] = int(g)
        #print('%s\t%d' % (q,fIndex[ind]))
  
  for i in range(nClusters):
    fcseqs  = []
    fcName = 'frag.' + str(i) + '.fas'
    
    for j in range(len(fSeqIDs)):
      if fIndex[j] == i:
        fcseqs.append(fSeqs[j])    
    
    if len(fcseqs) > 0:
      SeqIO.write(fcseqs, fcName, 'fasta')
      print('\t<%s> contains %d fragments' % (fcName,len(fcseqs)))
    
  #print(fSeqIDs)
  #print(fIndex) 
  
  noCluster = []
  
  ncIds = '' # gets the list of sequences not in any clusters
    
  for i in range(len(fSeqs)):
    if fIndex[i] == -1:
      noCluster.append(fSeqs[i])
      ncIds += fSeqs[i].id + '\t\t'
  
  if len(noCluster) > 0:
    SeqIO.write(noCluster,'frag.noClusters.fas','fasta')  
    print('\nThe following fragments did not match with any cluster')
    print('\n%s\n' % ncIds)

    if keepOrphans:
      print('But will be added to the final alignment')
      #print('\nUse the flag "-f 0" to exclude these from the final alignment')
      return True
    
    else: # user did not explicitely asked to add orphan fragments  
      print('**You can use "-f" flag to explicitely include fragments without clusters**\n')
      while True:
        choice = input('Do you want to add these sequences to the alignment (y/n)? ')
        if choice.strip() == 'y' or choice.strip() == 'Y':
          return True
        elif choice.strip() == 'n' or choice.strip() == 'N':
          return False
        else:
          print('\nIncorrect choice. Please enter "y(Y)" or "n(N)"')
          continue
    
#************************************************************************

#************************************************************************
def addFragmentsToClusters(nClusters,thread,cDir,tName,zName):
  '''
    Fragments are added to their corresponding cluster alignments
    MAFFT's "--addfragments" is used for adding fragments
  '''    

  log = 'fragAlign.log'
  lh = open(log,'w')
  lh.close()
  
  for i in range(nClusters):
    fName = 'frag.' + str(i) + '.fas'
    aName = 'grp.' + str(i) + '.aln'
    oName = 'cls.' + str(i) + '.aln'
    if os.path.exists(aName) and os.stat(aName).st_size > 0:
      if os.path.exists(fName) and os.stat(fName).st_size > 0:
        print('\nAdding fragments from <%s> to <%s>' % (fName,aName)) 
        addFragments(fName,aName,oName,thread,log,cDir,tName,zName)
      else:
        try:
          shutil.copy(aName,oName)
        except OSError as e:
          print(e)
          cZip(cDir,tName,zName)

  lh.close()
#************************************************************************

#***********************************************************************
def mergeClusters(nClusters,outFile,keepOrphans,thread,mIterM,cDir,tName,zName,lEx):
  '''
    - Merge clusters into one large alignment
    - adds the fragments that were not assigned any cluster if chosen by the user 
  '''  
  
  print('\nMerging all cluster alignments together')
  lh = open('merge.log','w')
  if nClusters == 1: 
    if not keepOrphans: 
      try:
        shutil.copy('cls.0.aln',outFile)
        return
      except OSError as e:
        print(e)
        cZip(cDir,tName,zName)
    else:
      if os.path.exists('frag.noClusters.fas') and os.stat('frag.noClusters.fas').st_size > 0:
        addFragments('frag.noClusters.fas','cls.0.aln',outFile,thread,'alignOrphan.log',cDir,tName,zName)
        return
      else:
        try:
          shutil.copy('cls.0.aln',outFile)
          return
        except OSError as e:
          print(e)
          cZip(cDir,tName,zName)
      
  seqCount = 1
  catText = 'cat '
  mTab = ''
  
  #*********
  oSeqs = [] # contains orphan long sequences and fragments
  mFlag = False
  fFlag = False
  
  for i in range(nClusters):
    if i in lEx:
      continue
    
    cName = 'cls.' + str(i) + '.aln' # name of the cluster alignment file
    seqs = list(SeqIO.parse(cName,'fasta')) # reading the alignment
    if len(seqs) > 1: # if more than one sequence in the alignment
      catText += cName + ' ' # add file name to the cat string 
      for k in range(len(seqs)): # write numbers of sequences in line
        mTab += str(seqCount) + ' '
        seqCount += 1
      mTab += '\n'
      mFlag = True
    else:
      oSeqs.append(seqs[0]) # if only one sequence in the alignment, write to orphan file

  if keepOrphans and os.path.exists('frag.noClusters.fas'): # if there exists orphan fragments
    feqs = list(SeqIO.parse('frag.noClusters.fas','fasta'))
    for seq in feqs:
      oSeqs.append(seq) # add fragments to orphan seqs file
      #print('test')
  
  if len(oSeqs) > 0: # if at least one orphan sequence present
    SeqIO.write(oSeqs,'orphanSeqs.fasta','fasta') # write orphan long and fragments in file
    catText += 'orphanSeqs.fasta '
  
  catText += '> merge'

  
    
  #*********
  if mFlag: # at least one alignment contains more than one sequences
    try:
      subprocess.check_call(catText,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      print(e)
      cZip(cDir,tName,zName)

    fh = open('subMSAtable','w')
    fh.write(mTab)
    fh.close()
  
    cl = 'mafft --preservecase --thread %d --localpair --maxiterate %d --merge subMSAtable merge > out.aln' % (thread,mIterM) # uses L-INS-i
    #cl = 'mafft --preservecase --thread %d --maxiterate %d --merge subMSAtable merge > out.aln' % (thread,mIterM) # uses FFT-NS-i
  
  elif len(oSeqs) > 1:
    cl = 'mafft --preservecase --thread %d --localpair --maxiterate %d orphanSeqs.fasta > out.aln' % (thread,mIterM) # if no alignment
    fFlag = True
  
  elif len(oSeqs) == 1:
    print('\nAligning only one sequence.')
    try:  
      shutil.copyfile('orphanSeqs.fasta',outFile)  
    except OSError as e:
      print(e)
      cZip(cDir,tName,zName)
    
  if mFlag or fFlag:
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      print(e)
      cZip(cDir,tName,zName)
  
    if os.path.exists('out.aln') and os.stat('out.aln').st_size > 0:
      try:
        shutil.copy('out.aln',outFile)
      except OSError as e:
        print(e)
        cZip(cDir,tName,zName)
  
  lh.close()
  
#************************************************************************
    
#************************************************************************  
def getArguments():
  '''
    Parses all the command line arguments from the user
  
  '''
  parser = argparse.ArgumentParser(description="Pipelign: creates multiple sequence alignment from FASTA formatted sequence file", formatter_class=argparse.RawTextHelpFormatter) #formatter_class=argparse.RawDescriptionHelpFormatter)  
  parser.add_argument('-i', '--input', required=True, help="Input sequence file")
  parser.add_argument('-o', '--output', required=True, help="Output alignment file")
  parser.add_argument('-t', '--thr', type=lengthThreshold, help="Length threshold for full sequences", default=0.5)
  parser.add_argument('-c', '--code', type=int, help="Genetic code for translation",default=1, choices=[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25])
  parser.add_argument('-a', '--alphabet', required=True, help='Input sequences can be DNA/Protein', choices=['dna','aa','rna'], default='dna')
  parser.add_argument('-f', '--keepOrphans', help='Add fragments without clusters', action="store_true")
  parser.add_argument('-z', '--mZip', help='Create zipped temporary files', action="store_true")
  parser.add_argument('-p', '--simPer', type=lengthThreshold, help="percent sequence similarity for clustering", default=0.8)
  parser.add_argument('-q', '--thread', nargs='?', const=1, type=int, help="Number of CPU to use for multithreads", default=1)
  parser.add_argument('-s', '--mIterateLong', type=iterThreshold, help="Number of iterations to refine long alignments", default=1000)
  parser.add_argument('-m', '--mIterateMerge', type=iterThreshold, help="Number of iterations to refine merged alignment", default=100)
  parser.add_argument('-d', '--tempDirPath', required=False, help="Path for temporary directory",default=None)
  parser.add_argument('-l', '--longSeqsOnly', help='Only align long sequences', action="store_true")
  parser.add_argument('-n', '--stage', type=int,default=5, choices=[1,2,3,4],
    help=textwrap.dedent('''\
    1  Make clusters alignments of long sequences
    2  Build HMMs of long cluster alignments
    3  Add fragments to cluster alignments
    4  Make final alignment
  	'''))

  parser.add_argument('-x', '--excludeClusters', help='Exclude clusters from final alignment', action="store_true")  
  
  args = parser.parse_args()

  return args  
  
#****************************************************************************

#****************************************************************************

if __name__=="__main__":
  args = getArguments()
  
  #init(args)
  oFile = os.getcwd() + '/' + args.output
  mArgs = MyStruct(
      inFile = args.input,
      outFile = oFile,
      lenThr = args.thr,
      gCode = args.code,
      alphabet = args.alphabet,
      keepOrphans = args.keepOrphans,
      makeZip = args.mZip,
      simPer = args.simPer,
      thread = args.thread,
      mIterL = args.mIterateLong,
      mIterM = args.mIterateMerge,
      fragEmpty = 1,
      longName = 'long.fas',
      fragName = 'frag.fas',
      tempDirPath = args.tempDirPath,
      longSeqsOnly = args.longSeqsOnly,
      stage = args.stage,
      excludeClusters = args.excludeClusters)
  
  if mArgs.stage == 3 and mArgs.longSeqsOnly:
    print('\nOptions -n (--stage) and -l (--longSeqsOnly) cannot be used in the same command')
    print('Please use only one of them in a single run')
    sys.exit()
  
  cDir = os.getcwd() # save current working directory
  tName1 = 'in.fas'
  tName2 = 'input.fas'
  
  # check whether input file exists or exit
  
  if not os.path.exists(mArgs.inFile) or os.stat(mArgs.inFile).st_size == 0:
      msg = '\n\nError: input sequence file "%s"  could not be found' % mArgs.inFile
      msg += '\n\t**Please run Pipelign again with correct input file name**'
      msg += '\n\tPipelign is exiting\n'
      sys.exit(msg)
  
  # get name for the zipped temporary directory
  timeNow = time.strftime('%Y-%m-%d-%H%M%S')
  zName = 'pipelign.' + timeNow 

  # create temporary directory
  
  if mArgs.tempDirPath is None: # no path provided for temp directory
    try:
      tempDir = tempfile.TemporaryDirectory() # create temporary directory to hold intermediary files
      tName = tempDir.name
    except OSError as e:
      sys.exit('\nError: system could not create temporary directory. Please try again')
  
  else:
    if os.path.exists(mArgs.tempDirPath):
      tempDir = mArgs.tempDirPath + '/' + zName
      tName = tempDir
      try:
        os.mkdir(tempDir)
      except OSError as e:
        sys.exit('\nError: system could not create temporary directory. Please try again')
    else:
      sys.exit('\nError: Path for temporary directory does not exists. Please run again with correct path.')
      
  # copy input file inside the temporary directory
  tFileName = tName + '/' + tName1
  
  try:  
    shutil.copyfile(mArgs.inFile,tFileName)  
  except OSError as e:
      sys.exit('\nError: could not copy input file into temp directory. Please try again ')
  
  
  print('\nPipelign will align sequences in: <%s>' % mArgs.inFile)
  
  
  #change current working directory to the temp
  os.chdir(tName)
  deAlign(tName1, tName2) # removes any possible gaps from the sequence file
  
  '''
  # plot a sequence length distribution and save it as svg file
  #cl = 'Rscript %s/pipelignRscripts/lengthDistribution.R %s ../lengthDistribution.svg' % (cDir,tName2)
  cl = 'lengthDistribution.R %s ../lengthDistribution.svg' % tName2
  try:
    subprocess.check_call(cl,shell=True,stdout=None,stderr=None)
  except subprocess.CalledProcessError as e:
    #sys.exit(e)
    print(e)
    cZip(cDir,tName,zName)
  print('\t<lengthDistribution.svg> created')
  '''
  
  mArgs.fragEmpty = separateFullFragment(tName2, mArgs.lenThr, mArgs.longName, mArgs.fragName)
  
  runCDHIT(mArgs.longName, mArgs.alphabet, mArgs.simPer, mArgs.thread,cDir,tName,zName)
    
  numClusters, clsSize = makeClusters(mArgs.longName)
  
  if numClusters == 0:
    sys.exit('\nError: Something went wrong during clustering sequences. Please check the input file')
    
  addClusterNumberToReps('grp','clusterList.txt','clsReps.fas',clsSize)
  #'''
  makeClusterRepsAlignment('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,cDir,tName,zName)
  
  clsExclude = list()
  
  if numClusters > 2:
    makeIQTree('clsReps.aln',mArgs.thread,cDir,tName,zName)
    #drawTree('clsReps.aln.treefile')
    drawMidPointRootTree('clsReps.aln.treefile')
    print('All %d cluster(s) will be added to the final alignment' % numClusters)
    
    if mArgs.excludeClusters:
      inChoice = input('Clusters you want to exclude: ')
      clsL = inChoice.split()
      for c in clsL:
        if 0 <= int(c) <= numClusters:
          clsExclude.append(int(c))
          print('\tCluster %d will not be added to the final alignment' % int(c))
        else:
          print('\t%d is not a valid cluster number' % int(c))
  
  else:
    print('\nNumber of cluster representative(s) is %d. Phylogenetic tree can not be built' % numClusters)
    
  #print('\n\nNumber of clusters %d' % numClusters)
  alnFullSequenceClusters(numClusters, mArgs.thread,mArgs.mIterL,cDir,tName,zName)

  if mArgs.stage == 1: # stops after producing clusters of long sequences
    print('\nAlignments of long sequence clusters created')
    cZip(cDir,tName,zName)

  if mArgs.stage == 2:
    lFile = 'hmmer.log'
    oFile = 'hmm.out'
    makeHMMdb(numClusters,cDir,tName,zName,mArgs.thread,lFile,mArgs.alphabet)
    print('\nHMM files are created.')
    cZip(cDir,tName,zName)
      
  if not mArgs.longSeqsOnly: # add fragments to final alignment
    if not mArgs.fragEmpty:
      lFile = 'hmmer.log'
      oFile = 'hmm.out'
      makeHMMdb(numClusters,cDir,tName,zName,mArgs.thread,lFile,mArgs.alphabet)
      searchHMMdb(lFile,mArgs.thread,mArgs.fragEmpty,mArgs.alphabet,oFile,cDir,tName,zName)
    
      #numClusters = 2
      mArgs.keepOrphans = parseHMMsearchResult(numClusters,mArgs.fragName,oFile,mArgs.keepOrphans)
  
      # next is add fragments to cluster alignments  
      addFragmentsToClusters(numClusters,mArgs.thread,cDir,tName,zName)
    
    '''
    #cl = 'Rscript %s/pipelignRscripts/freqsLongFrags.R clsReps.aln.treefile clusterList.txt ' % cDir
    cl = 'freqsLongFrags.R clsReps.aln.treefile clusterList.txt ' 
    cl += 'hmm.out ../stat.frequency.long.svg ../stat.frequency.fragments.svg'
    
    lh = open('stat.frequency.log','w')
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      #sys.exit(e)
      print(e)
      cZip(cDir,tName,zName)
    print('\t<stat.frequency.long.svg> created')
    print('\t<stat.frequency.fragments.svg> created')
    lh.close()
    '''
    if mArgs.fragEmpty: # no fragments present, writing cluster alignment files
      for i in range(numClusters):
        try:
          aName = 'grp.' + str(i) + '.aln'
          oName = 'cls.' + str(i) + '.aln'
          
          shutil.copy(aName,oName)
        except OSError as e:
          print(e)
          cZip(cDir,tName,zName)
    
    if mArgs.stage == 3:
      print('\nFragments are added to the cluster alignments')
      cZip(cDir,tName,zName)
  
  else: # only align long sequences
    print("\n**Pipelign is running with '-l' flag. So only long sequences will be added to the final alignment**")
    for i in range(numClusters):
      try:
        aName = 'grp.' + str(i) + '.aln'
        oName = 'cls.' + str(i) + '.aln'

        shutil.copy(aName,oName)
      except OSError as e:
        print(e)
        cZip(cDir,tName,zName) 

    '''
    #cl = 'Rscript %s/pipelignRscripts/freqsLong.R clsReps.aln.treefile clusterList.txt ' % cDir
    cl = 'freqsLong.R clsReps.aln.treefile clusterList.txt ' 
    cl += '../stat.frequency.long.svg'
    
    lh = open('stat.frequency.log','w')
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      #sys.exit(e)
      print(e)
      cZip(cDir,tName,zName)
    print('\t<stat.frequency.long.svg> created')
    lh.close()
    '''
 
  mergeClusters(numClusters,mArgs.outFile,mArgs.keepOrphans,mArgs.thread,mArgs.mIterM,cDir,tName,zName,clsExclude)

  print('\nThe alignment is written in <%s>\n' % mArgs.outFile)
  
  if mArgs.makeZip:
    cZip(cDir,tName,zName)
  #'''
#************************************************************  

