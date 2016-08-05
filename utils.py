import sys, os, shutil, subprocess, argparse, textwrap

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
