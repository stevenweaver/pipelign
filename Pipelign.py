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


class Struct:
  '''
   - this class creates a structure of the parameter values needed for different functions.
   - parameters can be accessed as structName.itemName e.g. mArgs.gCode 
  '''
  def __init__(self, **kwargs):
    for k, v in kwargs.items():
      setattr(self,k,v)
      
class MyStruct(Struct):
  pass

#********************************************************************

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

def cZip(cDir,tName):
  '''
  creates a zip file of the temporary directory
  '''  
  os.chdir(cDir)
  
  zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S') 
  shutil.make_archive(zName,'zip',tName)
  print('\nArchive for all temporary files created in %s.zip\n' % zName)
  sys.exit()

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
def runCDHIT(longName, alphabet, per, thread):
  '''
    CD-HIT is used to group similar sequences together in clusters for alignment
    
  '''
   
  seqs = list(SeqIO.parse(longName,'fasta'))
  
  if len(seqs) < 2:
    return
  
  print('\nRunning CD-HIT/CD-HIT-EST to group long sequences into clusters based on sequence similarity')
  
  lh = open('cdhit.log','w') # create log file for cdhit
  
  if alphabet == 'dna':
    cl = 'cd-hit-est -c %f -n 5 -i %s -o grp -d 0 -T %d' % (per,longName,thread)
    #print(cl)
    
    try:
      subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
      sys.exit(e)
      
  elif alphabet == 'protein':
    cl = 'cd-hit -c %f -n 5 -i %s -o grp -d 0 -T %d' % (per,longName,thread)
    #print(cl)
    
    try:
      subprocess.check_call(cl, shell=True, stdout=lh, stderr=lh)
    except subprocess.CalledProcessError as e:
      sys.exit(e)
  
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
def alnFullSequenceClusters(nClusters, thread):
  '''
    Full sequences in each clusters will be aligned using L-INS-i/clustalo
  
  '''
  lh = open('clusterAlign.log','w') # open log file for cluster alignment
  
  print('\nAligning clusters')
  for i in range(nClusters):
    cName = 'grp.' + str(i) + '.fas'
    aName = 'grp.' + str(i) + '.aln'
    seqs = list(SeqIO.parse(cName,'fasta'))
    
    
    
    if len(seqs) > 1:
      print('\tAligning cluster %d of %d sequences' % (i+1,len(seqs)))
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
def makeHMMdb(nClusters,cDir,tName,thread,lFile):
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
    
    cl = 'hmmbuild --cpu %d %s %s' % (thread, hName, aName)
    try:
      subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
    except subprocess.CalledProcessError as e:
      #sys.exit(e)
      print(e)
      cZip(cDir,tName)
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
def searchHMMdb(lFile, thread, fragEmpty,alpha,res):
  '''
    HMM database is searched with the fragments to assign them a cliuster for alignment
  '''  
  lh = open(lFile,'a')
  
  if not fragEmpty:
    print('\nFragments will be searched against the HMM database to assign clusters')
    print('\tResults of the database search are being written on <%s>' % res)
    
    if alpha == 'dna':
      cl = 'nhmmscan --cpu %d --tblout %s --noali grp.hmm frag.fas' % (thread,res)
      try:
        subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
      except subprocess.CalledProcessError as e:
        sys.exit(e)
      
    elif alpha == 'protein':
      cl = 'hmmscan --cpu %d --tblout %s --noali grp.hmm frag.fas' % (thread,res)
      try:
        subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
      except subprocess.CalledProcessError as e:
        sys.exit(e)
  
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
      print('\nUse the flag "-f 0" to exclude these from the final alignment')
    
    else: # user did not explicitely asked to add orphan fragments  
      while True:
        choice = input('Do you want to add these sequences to the alignment (y/n)? ')
        if choice.strip() == 'y' or choice.strip() == 'Y':
          return 1
        elif choice.strip() == 'n' or choice.strip() == 'N':
          return 0
        else:
          print('\nIncorrect choice. Please enter "y(Y)" or "n(N)"')
          continue
    
#************************************************************************

#************************************************************************
def addFragmentsToClusters(nClusters, thread):
  '''
    Fragments are added to their corresponding cluster alignments
    MAFFT's "--addfragments" is used for adding fragments
  '''    
  #if not fragEmpty:
    #print('\nAdding fragments to their corresponding cluster alignments')
  
  lh = open('fragAlign.log','w')
  
  for i in range(nClusters):
    fName = 'frag.' + str(i) + '.fas'
    aName = 'grp.' + str(i) + '.aln'
    oName = 'cls.' + str(i) + '.aln'
    if os.path.exists(aName) and os.stat(aName).st_size > 0:
      if os.path.exists(fName) and os.stat(fName).st_size > 0:
        print('\nAdding fragments from <%s> to <%s>' % (fName,aName)) 
        cl = 'mafft --preservecase --thread %d --addfragments %s %s > %s' % (thread, fName, aName, oName)
        
        try:
          subprocess.check_call(cl,shell=True,stdout=lh,stderr=lh)
        except subprocess.CalledProcessError as e:
          sys.exit(e)
      else:
        shutil.copy(aName,oName)

  lh.close()
#************************************************************************

def mergeClusters(nClusters,outFile,addNClusters,thread):
  '''
    - Merge clusters into one large alignment
    - adds the fragments that were not assigned any cluster if chosen by the user 
  '''  
  
  print('\nMerging all cluster alignments together')
  
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
  parser.add_argument('-a', '--alphabet', help='Input sequences can be DNA/Protein', choices=['dna','protein'], default='dna')
  parser.add_argument('-f', '--keepFrag', nargs = '?', const = 1, type = int, help='Add fragments without clusters', default=1)
  parser.add_argument('-z', '--mZip', nargs = '?', const = 0, type = int, help='Create zipped temporary files', default=0)
  parser.add_argument('-p', '--simPer', nargs='?', const=0.8, type=float, help="percent sequence similarity for clustering", default=0.8)
  parser.add_argument('-q', '--thread', nargs='?', const=1, type=int, help="Number of CPU to use for multithreads", default=1)
  
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
      fragEmpty = 1,
      longName = 'long.fas',
      fragName = 'frag.fas')
  
  
  cDir = os.getcwd() # save current working directory
  tName1 = 'in.fas'
  tName2 = 'input.fas'
  
  #'''
  # create temporary directory
  try:
    tempDir = tempfile.TemporaryDirectory() # create temporary directory to hold intermediary files
    tFileName = tempDir.name + '/' + tName1 
    shutil.copyfile(mArgs.inFile,tFileName) # copy input file to temporary directory 
  
  except OSError:
    msg = '\n\nThe sequence file "%s"  cannot be found.' % mArgs.inFile
    msg += '\n**Please run Pipelign again with correct input file name**'
    msg += '\nPipelign is exiting.\n'
    sys.exit(msg)
  
  print('\nPipelign will align sequences in: <%s>' % mArgs.inFile)
  tName = tempDir.name
  
  #change current working directory to the temp
  os.chdir(tName)
  deAlign(tName1, tName2) # removes any possible gaps from the sequence file
  mArgs.fragEmpty = separateFullFragment(tName2, mArgs.lenThr, mArgs.longName, mArgs.fragName)
  
  runCDHIT(mArgs.longName, mArgs.alphabet, mArgs.simPer, mArgs.thread)
  
  numClusters = makeClusters(mArgs.longName)
  
  #print('\n\nNumber of clusters %d' % numClusters)
  alnFullSequenceClusters(numClusters, mArgs.thread)
  
  if not mArgs.fragEmpty:
    lFile = 'hmmer.log'
    oFile = 'hmm.out'
    makeHMMdb(numClusters,cDir,tName,mArgs.thread,lFile)
    searchHMMdb(lFile,mArgs.thread,mArgs.fragEmpty,mArgs.alphabet,oFile)
    
    #numClusters = 2
    addNC = parseHMMsearchResult(numClusters,mArgs.fragName,oFile,mArgs.keepFrag)
    # next is add fragments to cluster alignments  
    addFragmentsToClusters(numClusters,mArgs.thread)
  
  mergeClusters(numClusters,mArgs.outFile,addNC,mArgs.thread)
  print('\nThe alignment is written in %s' % mArgs.outFile)
  
  
  os.chdir(cDir)
  
  # create zipped file for temporary directory if -z 1 is used
  if args.mZip:
    zName = 'pipelign.' + time.strftime('%Y-%m-%d-%H%M%S') 
    shutil.make_archive(zName,'zip',tempDir.name)
    print('\nArchive for all temporary files created in %s.zip\n' % zName)

#************************************************************  

