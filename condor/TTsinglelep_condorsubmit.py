#!/usr/bin/python

import os
import re
import sys
import fileinput
import getopt
import subprocess
import socket
import datetime
files_per_job = 10

#Absolute path that precedes '/store...'
# The script checks whether it is an absolute path or not

if (socket.gethostname().find('fnal')>=0):
    brux=bool(False)
    print "Submitting jobs at FNAL"
    sePath='\'root://cmseos.fnal.gov/'         # stored on LPC
#    sePath='\'root://cmsxrootd.fnal.gov/'      # stored elsewhere
    storePath=''
    setupString='source \/cvmfs\/cms.cern.ch\/cmsset_default.csh'
    outPath='root://cmseos.fnal.gov/'          # output on LPC

else:
    print "Script not done for ", socket.gethostname()
    exit()

#Configuration options parsed from arguments
#Switching to getopt for compatibility with older python
try:
    opts, args = getopt.getopt(sys.argv[1:], "", ["useMC=", "sample=","dataType=","fileList=", "outDir=","shift=","submit=","json=","inputTar="]) #"jec=","jer=","btag=",
except getopt.GetoptError as err:
    print str(err)
    sys.exit(1)

useMC    = bool(False)
prefix   = str('None')
dataType = str('None')
fileList = str('None')
outDir   = str('None')
shift    = str('None')
json     = str('None')
submit   = bool(False)
#jec      = str('False')
#jer      = str('False')
#btag     = str('False')
changeJEC= bool(False)
tarfile  = str('None')

for o, a in opts:
    print o, a
    if o == "--useMC":
        checkUseMC = a in ['True', 'False']
        if not checkUseMC:
            print '--useMC must be True or False!'
            print 'Setting useMC to False'
        else:
            useMC = a
            print 'useMC = ' + useMC
    elif o == "--sample":   prefix   = str(a)
    elif o == "--dataType": dataType = str(a)
    elif o == "--fileList": fileList = str(a)
    elif o == "--outDir":   outDir   = str(a)
    elif o == "--shift":    shift    = str(a)
    elif o == "--inputTar": tarfile = str(a).split('.')[0]
    elif o == "--json":     json   = str(a)
#    elif o == "--jec":      jec   = str(a)
#    elif o == "--jer":      jer   = str(a)
#    elif o == "--btag":     btag   = str(a)
    elif o == "--submit":
        if a == 'True':     submit = True

if dataType == 'MuEl': dataType = 'ElMu'
checkDataType = dataType in ['ElEl', 'ElMu', 'MuMu']

if (not checkDataType and not useMC):
    print '--dataType for real data must be ElEl, ElMu or MuMu!'
    sys.exit(1)

relBase = os.environ['CMSSW_BASE']

files = fileList 

if (not os.path.isfile(files)):
    print 'File with input root files '+ fileList +' not found. Please give absolute path'
    sys.exit(1)

outputdir = os.path.abspath('.')
if str(outDir) != 'None': outputdir = outDir
else:                     print 'Output dir not specified. Setting it to current directory.'

if submit: print 'Will submit jobs to Condor'
else:      print 'Will not submit jobs to Condor'

dir = outputdir+'/'+shift+'/'+prefix
xrddir = outPath+dir[10:]

if os.path.exists(dir):
    print 'Output directory already exists!'
    sys.exit(1)
else: 
    if dir[:4] == '/eos':
        print 'Making outputDir via eos command'
        os.system('eos root://cmseos.fnal.gov/ mkdir -p '+dir[10:])
    else:
        os.system('mkdir -p '+dir)
    
def get_input(num, list):
    result = ''
    file_list = open(files)
    file_count = 0
    for line in file_list:
        if line.find('root')>0:
            file_count=file_count+1
            if file_count>(num-1) and file_count<(num+files_per_job):
                result=result+ sePath
		if (line.strip()[:6]=='/store'):
		    result=result+ storePath
                result=result+line.strip()+'\',\n'
    file_list.close()
    return result

j = 1
nfiles = 1

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
time='%i_%i_%i'%(cTime.hour,cTime.minute,cTime.second)
tempdir = '/uscms_data/d3/jmanagan/'+outputdir.split('/')[-1]+'_logs/'+shift+'/'+prefix
os.makedirs(tempdir)

py_file = open(tempdir+"/"+prefix+"_"+str(j)+".py","w")

print 'CONDOR work dir: '+tempdir

count = 0
file_list = open(files)
for line in file_list:
    if line.find('.root')>0:
        count = count + 1

file_list.close()

os.system('sed -e \'s/SETUP/'+setupString+'/g\' template.sh > '+tempdir+'/'+prefix+'.sh')

while ( nfiles <= count ):    

    # ADD YOUR CONFIG FILE HERE!!
    py_templ_file = open(relBase+'/src/LJMet/Com/condor_1lep_120715/ljmet_cfg.py')
    
    py_file = open(tempdir+'/'+prefix+'_'+str(j)+'.py','w')

    #Avoid calling the get_input function once per line in template
    singleFileList = get_input(nfiles, files)
    
    for line in py_templ_file:
        line=line.replace('CONDOR_ISMC',     useMC)
        line=line.replace('CONDOR_DATATYPE', dataType)
        line=line.replace('CONDOR_JSON',     json)
        line=line.replace('CONDOR_FILELIST', singleFileList)
        line=line.replace('CONDOR_OUTFILE',  prefix+'_'+str(j))
        if 'JECup' in shift  and 'JECup' in line: line=line.replace('False',  'True')
        if 'JECdown' in shift and 'JECdown' in line: line=line.replace('False',  'True')
        if 'JERup' in shift and 'JERup' in line: line=line.replace('False',  'True')
        if 'JERdown' in shift and 'JERdown' in line: line=line.replace('False',  'True')
        if 'BTAGup' in shift and 'BTagUncertUp' in line: line=line.replace('False',  'True')
        if 'BTAGdown' in shift and 'BTagUncertDown' in line: line=line.replace('False',  'True')
        py_file.write(line)
    py_file.close()

    jdl_templ_file = open(relBase+'/src/LJMet/Com/condor_1lep_120715/template.jdl')
    jdl_file       = open(tempdir+'/'+prefix+'_'+str(j)+'.jdl','w')
    
    for line in jdl_templ_file:
        line=line.replace('PREFIX',     prefix)
        line=line.replace('JOBID',      str(j))
        line=line.replace('OUTPUT_DIR', xrddir)
        line=line.replace('INPUTTAR',   tarfile)
        line=line.replace('TAR_FILE',   tarfile.split('/')[-1])
        jdl_file.write(line)
        
    jdl_file.close()

    j = j + 1
    nfiles = nfiles + files_per_job
    py_templ_file.close()

njobs = j - 1

if (submit):
    savedPath = os.getcwd()
#    os.system('xrdcp -R '+tempdir+' '+xrddir+'/logfiles/')
#    os.chdir(dir+'/logfiles/')
    os.chdir(tempdir)
    print njobs
    proc = subprocess.Popen(['condor_submit', prefix+'_1.jdl'], stdout=subprocess.PIPE) 
    (out, err) = proc.communicate()
    print "condor output:", out

    for k in range(2,njobs+1):
        os.system('condor_submit '+prefix+'_'+str(k)+'.jdl')

    os.chdir(savedPath)

#os.system('rm -r '+tempdir)

