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


inFile = '' # name of the input sequence file
outFile = '' # name of the output alignment file
thr = 0.0 # length thereshold to separate full sequences and fragments
genCode = 1 # genetic code to use
fragEmpty = 1 # flag to check whether there is any fragmemts in the input 
alpha = '' # sequence is DNA/Protein
addNClusters = 0  # flag for adding fragments without clusters
frag = 1 # ask users whether to keep fragments without clusters
zzip = 0 # writes zipped temporary files if zzip = 1
per = 0.8
thread = 1

#********************************************************************

#******************************************************************
def init(args):
  '''
    This function initializes names for:
      - input file
      - outfile file
    and values for:
      - minimum length for full length sequences
      - genetic code for transalation
  '''
  global inFile, outFile, thr, genCode, alpha, frag, zzip, per, thread
  
  inFile = args.input
  outFile = os.getcwd() + '/' + args.output
  thr = args.thr
  genCode = args.code
  alpha = args.alpha
  frag = args.frag
  zzip = args.zip
  per = args.per
  thread = args.thread
  
  print('\nInput sequence file for alignment is <%s>' % inFile)
  
#************************************************************************
  
#************************************************************************

def separateFullFragment():
  '''
    Reads in the input sequence file.
      - finds length of the longest sequence
      - find minimum length required for full length sequences
      - writes full length sequences and fragments into two separate files 
  '''
  global fragEmpty
  
  print('\nInput sequence file is being read')
  seqs = list(SeqIO.parse('input.fas','fasta'))
  
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
      
  if len(frag) > 0:
    SeqIO.write(frag, 'frag.fas','fasta')
    fragEmpty = 0
  
  SeqIO.write(full,'full.fas','fasta') 
#************************************************************************

#************************************************************************
def runCDHIT():
  '''
    CD-HIT is used to group similar sequences together in clusters for alignment
    
  '''
   
  seqs = list(SeqIO.parse('full.fas','fasta'))
  
  if len(seqs) < 2:
    return
  
  print('\nrunning CD-HIT/CD-HIT-EST to group input sequences into clusters based on sequence similarity')
  
  lh = open('cdhit.log','w') # create log file for cdhit
  
  if alpha == 'dna':
    cl = 'cd-hit-est -c %f -n 5 -i full.fas -o grp -d 0 -T %d' % (per,thread)
    #print(cl)
    
    try:
      subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
      sys.exit(e)
      
  elif alpha == 'protein':
    cl = 'cd-hit -c %f -n 5 -i full.fas -o grp -d 0 -T %d' % (per, thread)
    #print(cl)
    
    try:
      subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
      sys.exit(e)
  
  lh.close() # close log file   
#*************************************************************************

#*************************************************************************
def makeClusters():
  '''
    Separate files are created for each clusters
  '''
  
  print('\nWriting sequences to cluster files')
  
  fName = 'full.fas'
  
  seqs = list(SeqIO.parse(fName,'fasta')) # load all the full sequences
  
  if len(seqs) < 2:
    shutil.copy(fName,'grp.0.fas')
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
  
  return cls   
#***********************************************************************

#***********************************************************************
def alnFullSequenceClusters(nClusters):
  '''
    Full sequences in each clusters will be aligned using L-INS-i/clustalo
  
  '''
  lh = open('clusterAlign.log','w') # open log file for cluster alignment
  
  for i in range(nClusters):
    cName = 'grp.' + str(i) + '.fas'
    aName = 'grp.' + str(i) + '.aln'
    seqs = list(SeqIO.parse(cName,'fasta'))
    
    
    
    if len(seqs) > 1:
      cl = 'mafft --globalpair --thread %d --maxiterate 1000 --preservecase %s > %s' % (thread, cName, aName)
      #cl = 'clustalo -i %s -o %s' % (cName, aName)
      try:
        subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
      except subprocess.CalledProcessError as e:
        sys.exit(e)
    else:
      shutil.copyfile(cName,aName)

  lh.close() # close log file
#***********************************************************************  

#***********************************************************************  
def makeHMMdb(nClusters):
  '''
    Create profile HMMs for each of the full length alignments 
    Create HMM_DB from the pHMMs
     
  '''
  if fragEmpty:
    return # if there is no fragment sequences

  
  print('\nCreating profile HMMs for each of the full length alignments')
  
  lh = open('hmmer.log','a') # create log file for hmmer
  
  for i in range(nClusters):
    aName = 'grp.' + str(i) + '.aln'
    hName = 'grp.' + str(i) + '.hmm'
    
    #*****testing for failure of hmmbuild
    # hmmbuild fails when alignment does not have any gaps
    # solution: add a '-' at the end of each sequences ????
    aseqs = list(SeqIO.parse(aName,'fasta'))
    
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
    
    cl = 'hmmbuild --cpu %d %s %s' % (thread, hName, aName)
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      sys.exit(e)
    print('\n%s created' % hName)
  
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
    sys.exit(e)

  print('\nHMM database %s is created' % dName)
    
  lh.close() # close log file     
#***********************************************************************

#***********************************************************************  
def searchHMMdb():
  '''
    HMM database is searched with the fragments to assign them a cliuster for alignment
  '''  
  lh = open('hmmer.log','a')
  
  if not fragEmpty:
    print('\nFragments are to be searched')
    
    if alpha == 'dna':
      cl = 'nhmmscan --cpu %d --tblout hmm.out --noali grp.hmm frag.fas' % thread
      try:
        subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
      except subprocess.CalledProcessError as e:
        sys.exit(e)
      
    elif alpha == 'protein':
      cl = 'hmmscan --cpu %d --tblout hmm.out --noali grp.hmm frag.fas' % thread
      try:
        subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
      except subprocess.CalledProcessError as e:
        sys.exit(e)
  
  else:
    print('\nNo fragment sequence present')
    
  lh.close()
#**********************************************************************

#**********************************************************************    
def parseHMMsearchResult(nClusters):
  '''
    Reads in the result file from hmm.out file to determine cluster for fragments
  '''
  
  global addNClusters
  
  fSeqIDs = [] # IDs of the fragments
  
  fFlag = 0
  
  if fragEmpty:
    return # if there is no fragment sequences
        
  fSeqs = list(SeqIO.parse('frag.fas','fasta'))
  
  for seq in fSeqs:
    fSeqIDs.append(seq.id)
  
  fIndex = [-1 for x in range(len(fSeqIDs))] # to hold cluster number assigned for each fragments
  
  lines = [line.strip('\n') for line in open('hmm.out','r')]
  
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
    print('\n\nThe following Sequences did not match with any cluster')
    print('\n%s\n' % ncIds)
    print('\nThese sequences will not be added to the final alignment')
    
    #addNClusters = 1
    
    if frag == 1: # asks the user whether to add the fragments that did not have clusters
      while True:
        res = input('\nDo you want to add these sequences to the alignment (y/n)? ')
        if res == 'y' or res == 'Y':
          addNClusters = 1
          break
        elif res == 'n' or res == 'N':
          break
        else:
          print('\nIncorrect choice. Please enter "y/Y" or "n/N"')
          continue
    
#************************************************************************

#************************************************************************
def addFragmentsToClusters(nClusters):
  '''
    Fragments are added to their corresponding cluster alignments
    MAFFT's "--addfragments" is used for adding fragments
  '''    
  if not fragEmpty:
    print('\nAdding fragments to their corresponding cluster alignments')
  
  lh = open('fragAlign.log','w')
  
  for i in range(nClusters):
    fName = 'frag.' + str(i) + '.fas'
    aName = 'grp.' + str(i) + '.aln'
    oName = 'cls.' + str(i) + '.aln'
    if os.path.exists(aName) and os.stat(aName).st_size > 0:
      if os.path.exists(fName) and os.stat(fName).st_size > 0:
        #print('\n%s\t%s' % (aName,fName)) 
        cl = 'mafft --preservecase --thread %d --addfragments %s %s > %s' % (thread, fName, aName, oName)
        
        try:
          subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
        except subprocess.CalledProcessError as e:
          sys.exit(e)
      else:
        shutil.copy(aName,oName)

  lh.close()
#************************************************************************

def mergeClusters(nClusters):
  '''
    - Merge clusters into one large alignment
    - adds the fragments that were not assigned any cluster if chosen by the user 
  '''  
  
  print('\nMerging the clusters into a single alignment')
  
  if nClusters < 2 and addNClusters == 0: 
    shutil.copy('cls.0.aln',outFile)
    return
  
  seqCount = 1
  catText = 'cat '
  mTab = ''
  
  for i in range(nClusters):
    cName = 'cls.' + str(i) + '.aln'
    catText += cName + ' '
    seqs = list(SeqIO.parse(cName,'fasta'))
    for k in range(len(seqs)):
      mTab += str(seqCount) + ' '
      seqCount = seqCount + 1
  
    mTab += '\n'   
  #print(mTab)
  
  if addNClusters == 1 and os.path.exists('frag.noClusters.fas'):
    catText += 'frag.noClusters.fas '  
  
  catText += '> merge'
  #print(catText)
  
  lh = open('merge.log','w')
  try:
    subprocess.check_call(catText,shell=True,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    sys.exit(e)

  fh = open('subMSAtable','w')
  fh.write(mTab)
  fh.close()
  
  
  cl = 'mafft --preservecase --thread %d --localpair --maxiterate 100 --merge subMSAtable merge > out.aln' % thread
  
  try:
    subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
  except subprocess.CalledProcessError as e:
    sys.exit(e)
  
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
  parser.add_argument('-t', '--thr', nargs='?', const=0.5, type=float, help="Length threshold for full sequences", default=0.5)
  parser.add_argument('-c', '--code', nargs='?', const=1, type=int, help="Genetic code for translation",default=1)
  parser.add_argument('-a', '--alpha', help='Input sequences can be DNA/Protein', choices=['dna','protein'], default='dna')
  parser.add_argument('-f', '--frag', nargs = '?', const = 1, type = int, help='Add fragments without clusters', default=1)
  parser.add_argument('-z', '--zip', nargs = '?', const = 0, type = int, help='Create zipped temporary files', default=0)
  parser.add_argument('-p', '--per', nargs='?', const=0.8, type=float, help="percent sequence similarity for clustering", default=0.8)
  parser.add_argument('-q', '--thread', nargs='?', const=1, type=int, help="Number of CPU to use for multithreads", default=1)
  
  args = parser.parse_args()

  return args  
  
#****************************************************************************

#****************************************************************************

if __name__=="__main__":
  args = getArguments()
  
  init(args)
  
  cDir = os.getcwd() # save current working directory
  
  # create temporary directory
  try:
    tempDir = tempfile.TemporaryDirectory() # create temporary directory to hold intermediary files
    tFileName = tempDir.name + '/input.fas' 
    shutil.copyfile(inFile,tFileName) # copy input file to temporary directory 
  
  except OSError:
    msg = '\n\nThe sequence file "%s"  cannot be found.' % inFile
    msg += '\n**Please run Pipelign again with correct input file name**'
    msg += '\nPipelign is exiting.\n'
    sys.exit(msg)
  
  tName = tempDir.name
  
  #change current working directory to the temp
  os.chdir(tName)
  separateFullFragment()
  runCDHIT()
  numClusters = makeClusters()
  #print('\n\nNumber of clusters %d' % numClusters)
  alnFullSequenceClusters(numClusters)
  makeHMMdb(numClusters)
  searchHMMdb()
  #numClusters = 2
  parseHMMsearchResult(numClusters)
  # next is add fragments to cluster alignments
  addFragmentsToClusters(numClusters)
  mergeClusters(numClusters)
  print('\nThe alignment is written in %s' % outFile)
  
  os.chdir(cDir)
  
  # create zipped file for temporary directory if -z 1 is used
  if args.zip:
    zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S') 
    shutil.make_archive(zName,'zip',tempDir.name)
    print('\nArchive for all temporary files created in %s.zip\n' % zName)
#************************************************************  

