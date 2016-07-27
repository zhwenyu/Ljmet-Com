#!/usr/env/python
import os

datasets=[
"/DoubleMuon/Run2016B-PromptReco-v1/MINIAOD",
"/DoubleMuon/Run2016B-PromptReco-v2/MINIAOD",
"/DoubleMuon/Run2016C-PromptReco-v2/MINIAOD",
"/DoubleMuon/Run2016D-PromptReco-v2/MINIAOD",
"/DoubleEG/Run2016B-PromptReco-v1/MINIAOD",
"/DoubleEG/Run2016B-PromptReco-v2/MINIAOD",
"/DoubleEG/Run2016C-PromptReco-v2/MINIAOD",
"/DoubleEG/Run2016D-PromptReco-v2/MINIAOD",
"/MuonEG/Run2016B-PromptReco-v1/MINIAOD",
"/MuonEG/Run2016B-PromptReco-v2/MINIAOD",
"/MuonEG/Run2016C-PromptReco-v2/MINIAOD",
"/MuonEG/Run2016D-PromptReco-v2/MINIAOD",
]

files = []

for dataset in datasets:
    #make filename from datasetname, first remove leading slash
    fname = dataset[1:len(dataset)]
    #then remove '/MINIAODSIM'
    fname = fname.replace('/MINIAOD','')
    #change middle slash to underscor
    fname = fname.replace('/','_')
    #add .txt to end
    fname = fname+'.txt'
    files.append(fname)
    #dump files
    command = "python das_client.py --query=\"file dataset=%s\" --limit=20000 >& %s" % (dataset,fname)
    print command
    os.system(command)


for file in files:
    f = open(file,'r')
    lines = f.readlines()
    f.close()
    f = open(file,'w')
    for line in lines:
        if line.find('/store')!=-1:
            f.write(line)
    f.close()



    
