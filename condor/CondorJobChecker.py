#/usr/env/python

import os, argparse


files=[]
MCdirectories=os.listdir("/uscms_data/d3/clint/using_git/T53/ljmet/CMSSW_7_4_10_patch1/src/LJMet/Com/python/Samples_Spring15MC")
Data15Bdirectories=os.listdir("/uscms_data/d3/clint/using_git/T53/ljmet/CMSSW_7_4_10_patch1/src/LJMet/Com/python/Samples_Run2015B")

dirs=[]

for dir in MCdirectories:
    dirs.append(dir)

for dir in Data15Bdirectories:
    dirs.append(dir)

file='/uscms_data/d3/clint/using_git/T53/ljmet/CMSSW_7_4_10_patch1/src/LJMet/Com/python/Samples_Spring15MC/'

for dir in dirs:

    if dir.find("TTJets")>=0 and dir.find("50ns")>=0 and dir.find(".txt~")==-1:
        file+=dir

print file

fileIn=open(file,'r')
n=0
for f in fileIn:
    n+=1

outfiles=os.listdir(outDirs{sample})
