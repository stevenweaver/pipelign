'''
  Pipelign.py

  A python based program to align virus sequences.

  The program takes as input a single FASTA formatted file

  Returns a FASTA formatted alignment file

'''
#*********************************************************************

import sys, os, shutil, subprocess, argparse
from Bio import SeqIO
import tempfile
import time

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

#************************************************************************

def separateFullFragment(iFile, thr, longName, fragName):
  '''
    Reads in the input sequence file.
      - finds length of the longest sequence
      - find minimum length required for full length sequences
      - writes full length sequences and fragments into two separate files 
  '''
  #global fragEmpty
  
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
    fragEmpty = 0
    
  
  if fragEmpty:
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
    return
  
  print('\nRunning CD-HIT/CD-HIT-EST to group long sequences into clusters based on sequence similarity')
  
  lh = open('cdhit.log','w') # create log file for cdhit
  
  if alphabet == 'nt':
    cl = 'cd-hit-est -c %f -n 5 -i %s -o grp -d 0 -T %d' % (per,longName,thread)
    #print(cl)
    
    try:
      subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
      print(e)
      cZip(cDir,tName,zName)
      
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
  
  if len(seqs) < 2:
    shutil.copy(longName,'grp.0.fas')
    return 1
  
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
  
  return cls   
#***********************************************************************

#***********************************************************************

def addClusterNumberToReps(repName,lstFile,outFile):
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
      
      seq.id = seq.id + '_' + cNum[ind]
      seq.name = seq.id
      seq.description = seq.id
    else:
      sys.exit('\nSequence %s does not have a cluster. Pipelign is exiting' % seq.id)
  
  SeqIO.write(seqs,outFile,'fasta')    

#***********************************************************************

#***********************************************************************

def addFragments(fName,aName,oName,thread,log,cDir,tName,zName):
  '''
    - add fragments to the alignment using MAFFT's --add
  '''

  lh = open(log,'a')
  
  cl = 'mafft --preservecase --thread %d --addfragments %s %s > %s' % (thread, fName, aName, oName)
        
  try:
    subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  lh.close()
#***********************************************************************

#***********************************************************************

def makeClusterRepsAlignment(repFile,outFile,thread,mIterL,cDir,tName,zName):
  '''
    - Creates a multiple sequence alignment from cluster reps with cluster numbers
    - uses G-INS-i
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
    
    
    
    if len(seqs) > 1:
      print('\tAligning cluster %d of %d sequences' % (i+1,len(seqs)))
      #ginsi(cName,aName,thread,mIterL,log,cDir,tName,zName)
      linsi(cName,aName,thread,mIterL,log,cDir,tName,zName)
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
      aseqs[1].id = 'temp'

    if len(aseqs) > 1:
      sumGap = 0
      for seq in aseqs:
        sumGap = sumGap + seq.seq.count('-')
      if sumGap == 0:
        for seq in aseqs:
          seq.seq = seq.seq + '-'
        SeqIO.write(aseqs,'temp.aln','fasta')  
        aName = 'temp.aln'
    
    #**************
    if alpha == 'nt':
      cl = 'hmmbuild --dna --cpu %d %s %s' % (thread, hName, aName)
    
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
    
    if alpha == 'nt':
      cl = 'nhmmscan --cpu %d --tblout %s --noali grp.hmm frag.fas' % (thread,res)
      try:
        subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
      except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
      
    elif alpha == 'aa':
      cl = 'hmmscan --cpu %d --tblout %s --noali grp.hmm frag.fas' % (thread,res)
      try:
        subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
      except subprocess.CalledProcessError as e:
        print(e)
        cZip(cDir,tName,zName)
  
  #else:
    #print('\nNo fragment sequence present')
    
  lh.close()
#**********************************************************************

#**********************************************************************    
def parseHMMsearchResult(nClusters,fragFile,res,keepFrag):
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
      print('\t<%s>' % fcName)
    
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

    if keepFrag:
      print('But will be added to the final alignment')
      #print('\nUse the flag "-f 0" to exclude these from the final alignment')
      return True
    
    else: # user did not explicitely asked to add orphan fragments  
      print('Please use "-f" flag to include fragments without clusters')
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
        shutil.copy(aName,oName)

  lh.close()
#************************************************************************

def mergeClusters(nClusters,outFile,addNClusters,thread,mIterM,cDir,tName,zName):
  '''
    - Merge clusters into one large alignment
    - adds the fragments that were not assigned any cluster if chosen by the user 
  '''  
  
  print('\nMerging all cluster alignments together')
  
  if nClusters == 1 and not addNClusters: 
    shutil.copy('cls.0.aln',outFile)
    return
  
  seqCount = 1
  catText = 'cat '
  mTab = ''
  
  '''
  for i in range(nClusters):
    cName = 'cls.' + str(i) + '.aln'
    catText += cName + ' '
    seqs = list(SeqIO.parse(cName,'fasta'))
    for k in range(len(seqs)):
      mTab += str(seqCount) + ' '
      seqCount = seqCount + 1
  
    mTab += '\n'   
  #print(mTab)
  
  if addNClusters and os.path.exists('frag.noClusters.fas'):
    catText += 'frag.noClusters.fas ' 
  
  catText += '> merge'
  #print(catText)
  '''
  #*********
  oSeqs = [] # contains orphan long sequences and fragments
  mFlag = False
  
  for i in range(nClusters):
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

  if addNClusters and os.path.exists('frag.noClusters.fas'): # if there exists orphan fragments
    feqs = list(SeqIO.parse('frag.noClusters.fas','fasta'))
    for seq in feqs:
      oSeqs.append(seq) # add fragments to orphan seqs file
      #print('test')
  
  if len(oSeqs) > 0: # if at least one orphan sequence present
    SeqIO.write(oSeqs,'orphanSeqs.fasta','fasta') # write orphan long and fragments in file
    catText += 'orphanSeqs.fasta ' 
  
  catText += '> merge'

  lh = open('merge.log','w')
    
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
  
  else:
    cl = 'mafft --preservecase --thread %d --localpair --maxiterate %d orphanSeqs.fasta > out.aln' % (thread,mIterM) # if no alignment
    
  try:
    subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    print(e)
    cZip(cDir,tName,zName)
  
  if os.path.exists('out.aln') and os.stat('out.aln').st_size > 0:
    shutil.copy('out.aln',outFile)
    #print(outFile)
  
  lh.close()
  
#************************************************************************
    
#************************************************************************  
def getArguments():
  '''
    Parses all the command line arguments from the user
  
  '''
  parser = argparse.ArgumentParser(description="Pipelign: creates multiple sequence alignment from FASTA formatted sequence file", formatter_class=argparse.RawDescriptionHelpFormatter)  
  parser.add_argument('-i', '--input', required=True, help="Input sequence file")
  parser.add_argument('-o', '--output', required=True, help="Output alignment file")
  parser.add_argument('-t', '--thr', type=lengthThreshold, help="Length threshold for full sequences", default=0.5)
  parser.add_argument('-c', '--code', type=int, help="Genetic code for translation",default=1, choices=[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25])
  parser.add_argument('-a', '--alphabet', required=True, help='Input sequences can be DNA/Protein', choices=['nt','aa'], default='nt')
  parser.add_argument('-f', '--keepFrag', help='Add fragments without clusters', action="store_true")
  parser.add_argument('-z', '--mZip', help='Create zipped temporary files', action="store_true")
  parser.add_argument('-p', '--simPer', nargs='?', const=0.8, type=float, help="percent sequence similarity for clustering", default=0.8)
  parser.add_argument('-q', '--thread', nargs='?', const=1, type=int, help="Number of CPU to use for multithreads", default=1)
  parser.add_argument('-s', '--mIterateLong', nargs='?', const=1000, type=int, help="Number of iterations to refine long alignments", default=1000)
  parser.add_argument('-m', '--mIterateMerge', nargs='?', const=100, type=int, help="Number of iterations to refine merged alignment", default=100)
  parser.add_argument('-d', '--tempDirPath', required=False, help="Path for temporary directory",default=None)
  
  
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
      keepFrag = args.keepFrag,
      makeZip = args.mZip,
      simPer = args.simPer,
      thread = args.thread,
      mIterL = args.mIterateLong,
      mIterM = args.mIterateMerge,
      fragEmpty = 1,
      longName = 'long.fas',
      fragName = 'frag.fas',
      tempDirPath = args.tempDirPath)
  
  
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
  zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S') 

  # create temporary directory
  
  if mArgs.tempDirPath is None: # no path provided for temp directory
    try:
      tempDir = tempfile.TemporaryDirectory() # create temporary directory to hold intermediary files
      tName = tempDir.name
    except OSError:
      sys.exit('\nError: system could not create temporary directory. Please try again')
  
  else:
    if os.path.exists(mArgs.tempDirPath):
      tempDir = mArgs.tempDirPath + '/' + zName
      tName = tempDir
      try:
        os.mkdir(tempDir)
      except OSError:
        sys.exit('\nError: system could not create temporary directory. Please try again')
    else:
      sys.exit('\nError: Path for temporary directory does not exists. Please run again with correct path.')
      
  # copy input file inside the temporary directory
  tFileName = tName + '/' + tName1
  
  try:  
    shutil.copyfile(mArgs.inFile,tFileName)  
  except OSError:
      sys.exit('\nError: could not copy input file into temp directory. Please try again ')
  
  
  print('\nPipelign will align sequences in: <%s>' % mArgs.inFile)
  
  
  #change current working directory to the temp
  os.chdir(tName)
  deAlign(tName1, tName2) # removes any possible gaps from the sequence file
  mArgs.fragEmpty = separateFullFragment(tName2, mArgs.lenThr, mArgs.longName, mArgs.fragName)
  
  runCDHIT(mArgs.longName, mArgs.alphabet, mArgs.simPer, mArgs.thread,cDir,tName,zName)
    
  numClusters = makeClusters(mArgs.longName)
  
  addClusterNumberToReps('grp','clusterList.txt','clsReps.fas')
  
  makeClusterRepsAlignment('clsReps.fas','clsReps.aln',mArgs.thread,mArgs.mIterL,cDir,tName,zName)
  
  if numClusters > 2:
    makeIQTree('clsReps.aln',mArgs.thread,cDir,tName,zName)
  
  else:
    print('\nNumber of cluster representative(s) is %d. Phylogenetic tree can not be built' % numClusters)
    
  #print('\n\nNumber of clusters %d' % numClusters)
  alnFullSequenceClusters(numClusters, mArgs.thread,mArgs.mIterL,cDir,tName,zName)
  
  if not mArgs.fragEmpty:
    lFile = 'hmmer.log'
    oFile = 'hmm.out'
    makeHMMdb(numClusters,cDir,tName,zName,mArgs.thread,lFile,mArgs.alphabet)
    searchHMMdb(lFile,mArgs.thread,mArgs.fragEmpty,mArgs.alphabet,oFile,cDir,tName,zName)
    
    #numClusters = 2
    addNC = parseHMMsearchResult(numClusters,mArgs.fragName,oFile,mArgs.keepFrag)
    # next is add fragments to cluster alignments  
    addFragmentsToClusters(numClusters,mArgs.thread,cDir,tName,zName)
  
  mergeClusters(numClusters,mArgs.outFile,addNC,mArgs.thread,mArgs.mIterM,cDir,tName,zName)
  print('\nThe alignment is written in %s' % mArgs.outFile)
  
  if mArgs.makeZip:
    cZip(cDir,tName,zName)
#************************************************************  

