import Pipelign
import task_cdhit,task_tree
import pytest

#def test_aln_full_sequence_clusters():
#    nClusters = 15
#    thread = 64
#    mIterL = 1000
#    cDir = './'
#    tName = './tests/rsrc/'
#    zName = 'pipelign.2016-07-28-172850'
#    Pipelign.alnFullSequenceClusters(nClusters, thread, mIterL, cDir, tName, zName)

#def test_cd_hit():
#    nClusters = 15
#    thread =32
#    mIterL = 1000
#    cDir = './'
#    tName = './tests/rsrc/'
#    zName = 'pipelign.2016-07-28-172850'
#    longName = './tests/rsrc/input.fas'
#    c = task_cdhit.runCDHIT(longName=longName, alphabet='dna', per=0.8, thread=thread, cDir=cDir, tName=tName, zName=zName).run()

#def test_make_clusters():
#    longName = './tests/rsrc/input.fas'
#    c = task_cdhit.makeClusters(longName=longName).run()

#def test_add_cluster_number_to_reps():
#    repName = 'grp'
#    lstFile = 'clusterList.txt'
#    outFile = 'clsReps.fas'
#    clsSize = 15
#    c = task_cdhit.addClusterNumberToReps(
#                    repName=repName,
#                    lstFile=lstFile,
#                    outFile=outFile,
#                    clsSize=clsSize).run()

#def test_add_cluster_number_to_reps():
#    repName = 'grp'
#    lstFile = 'clusterList.txt'
#    outFile = 'clsReps.fas'
#    clsSize = 15
#    c = task_cdhit.addClusterNumberToReps(
#                    repName=repName,
#                    lstFile=lstFile,
#                    outFile=outFile,
#                    clsSize=clsSize).run()

#def test_make_iq_tree():
#    repName = 'grp'
#    lstFile = 'clusterList.txt'
#    outFile = 'clsReps.fas'
#    clsSize = 15
#    c = task_cdhit.addClusterNumberToReps(
#                    repName=repName,
#                    lstFile=lstFile,
#                    outFile=outFile,
#                    clsSize=clsSize).run()

def test_draw_midpoint_root_tree():
    treeFile =
    c = task_tree.drawMidPointRootTree(treeFile=treeFile).run()

#def test_searchhmmdb():
#    nClusters = 15
#    thread = 64
#    mIterL = 1000
#    cDir = '/home/sweaver/programming/sfrost/'
#    tName = './rsrc/'
#    zName = 'pipelign.2016-07-28-172850'
#    Pipelign.searchHMMdb(lFile,thread,fragEmpty,alpha,res,cDir,tName,zName)
