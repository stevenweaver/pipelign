import Pipelign
import pytest

def test_aln_full_sequence_clusters():
    nClusters = 15
    thread = 64
    mIterL = 1000
    cDir = './'
    tName = './tests/rsrc/'
    zName = 'pipelign.2016-07-28-172850'
    Pipelign.alnFullSequenceClusters(nClusters, thread, mIterL, cDir, tName, zName)

#def test_searchhmmdb():
#    nClusters = 15
#    thread = 64
#    mIterL = 1000
#    cDir = '/home/sweaver/programming/sfrost/'
#    tName = './rsrc/'
#    zName = 'pipelign.2016-07-28-172850'
#    Pipelign.searchHMMdb(lFile,thread,fragEmpty,alpha,res,cDir,tName,zName)
