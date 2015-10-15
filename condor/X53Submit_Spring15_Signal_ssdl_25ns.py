#!/usr/bin/python

import os
import re
import fileinput

files_per_job = 1

rel_base = os.environ['CMSSW_BASE']

outdir = '/eos/uscms/store/user/lpctlbsm/clint/Spring15/25ns/Oct06/'

### What is the name of your FWLite Analyzer
FWLiteAnalyzer = 'ljmet'

### Which Systematics to do
DONOMINAL = 'True'
DOJES = 'False'
DOJER = 'False'
DOBTAG = 'False'
DOQCDMC = 'False'
DOTTBARSYS = 'False'

### JSON file to use
MYJSON = "''"

### Systematics flags
BTAGUNCERTUP = 'False'
BTAGUNCERTDOWN = 'False'
JECUNCERTUP = 'False'
JECUNCERTDOWN = 'False'
JERUNCERTUP = 'False'
JERUNCERTDOWN = 'False'


#################################################
### Names to give to your output root files
#################################################

prefix = []

if DONOMINAL=='True':
    prefix.extend([
            'X53X53To2L2Nu_LH_700',
            'X53X53To2L2Nu_RH_700',
            'X53X53To2L2Nu_LH_800',
            'X53X53To2L2Nu_RH_800',
            'X53X53To2L2Nu_LH_900',
            'X53X53To2L2Nu_RH_900',
            'X53X53To2L2Nu_LH_1000',
            'X53X53To2L2Nu_RH_1000',
            'X53X53To2L2Nu_LH_1100',
            'X53X53To2L2Nu_RH_1100',
            'X53X53To2L2Nu_LH_1200',
            'X53X53To2L2Nu_RH_1200',
            'X53X53To2L2Nu_LH_1300',
            'X53X53To2L2Nu_RH_1300',
            'X53X53To2L2Nu_LH_1400',
            'X53X53To2L2Nu_RH_1400',
            'X53X53To2L2Nu_LH_1500',
            'X53X53To2L2Nu_RH_1500',
            'X53X53To2L2Nu_LH_1600',
            'X53X53To2L2Nu_RH_1600',
    ])


###########################################
### Where to save your output root files to
###########################################
dir = []
for i in prefix:
    outdiri = outdir+i
    dir.extend([outdiri])

################################################
### Where is your list of root files to run over
################################################
list = [] 

listnom = [
'Samples_Spring15MC/X53X53To2L2Nu_M-700_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-700_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-800_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-900_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v2.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-900_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1100_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1100_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1200_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1300_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1300_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1500_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1500_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
'Samples_Spring15MC/X53X53To2L2Nu_M-1600_RH_TuneCUETP8M1_13TeV-madgraph-pythia8_RunIISpring15DR74-Asympt25nsRaw_MCRUN2_74_V9-v1.txt',
    ]

if DONOMINAL=='True':
    list.extend(listnom)
   

for i in range(len(prefix)):
    print i,' ',prefix[i],' ',dir[i],' ',list[i]

### Write the files you wish to run over for each job    
def get_input(num, list):
    result = '' 
    file_list = open(rel_base+"/src/LJMet/Com/python/"+list)
    file_count = 0
    for line in file_list:
        if line.find('root')>0:
            file_count=file_count+1
            if file_count>(num-1) and file_count<(num+files_per_job):
                f_name=line.split('.root')[0]
                f_name=f_name+'.root'
                #result=result+'                 \'' + f_name.group(1)+'\',\n'
                #result=result+'                 \'dcap:///pnfs/cms/WAX/11' + f_name.group(1)+'\',\n'
                result=result+'                 \'root://cmsxrootd.fnal.gov/' + f_name +'\',\n'
    file_list.close()
    #result = result + '                 )\n'
    return result


print str(files_per_job)+' files per job...'

for i in range(len(prefix)):

    j = 1
    nfiles = 1

        
    FLAGTAG = 'TriggerResults::PAT'
    
    print 'CONDOR work dir: '+dir[i]
    #os.system('rm -rf '+dir[i])
    os.system('mkdir -p '+dir[i])

    file_list = open(rel_base+"/src/LJMet/Com/python/"+list[i])
    count = 0
    for line in file_list:
        if line.find('root')>0:
            count = count + 1
    file_list.close()
    #count = count - 1

    print 'File prefix: '+prefix[i]
    print 'Number of input files: '+str(count)

    while ( nfiles <= count ):    

        py_templ_file = open(rel_base+"/src/LJMet/Com/condor/Dilepton_Spring15MC_25ns_python.templ")
        condor_templ_file = open(rel_base+"/src/LJMet/Com/condor/X53condor.templ")
        csh_templ_file    = open(rel_base+"/src/LJMet/Com/condor/X53csh.templ")

        py_file = open(dir[i]+"/"+prefix[i]+"_"+str(j)+".py","w")
        for line in py_templ_file:
            line=line.replace('DIRECTORY',dir[i])
            line=line.replace('PREFIX',prefix[i])
            line=line.replace('JOBID',str(j))
            line=line.replace('INFILES',get_input(nfiles, list[i]))
            line=line.replace('BTAGUNCERTUP',BTAGUNCERTUP)
            line=line.replace('BTAGUNCERTDOWN',BTAGUNCERTDOWN)
            line=line.replace('JECUNCERTUP',JECUNCERTUP)
            line=line.replace('JECUNCERTDOWN',JECUNCERTDOWN)
            line=line.replace('JERUNCERTUP',JERUNCERTUP)
            line=line.replace('JERUNCERTDOWN',JERUNCERTDOWN)
            line=line.replace('EVENTSTOPROCESS',str(-1))
            line=line.replace('FLAGTAG',FLAGTAG)
            line=line.replace('JSONFILE',MYJSON)
            py_file.write(line)
        py_file.close()

        condor_file = open(dir[i]+"/"+prefix[i]+"_"+str(j)+".condor","w")
        for line in condor_templ_file:
            line=line.replace('DIRECTORY',dir[i])
            line=line.replace('PREFIX',prefix[i])
            line=line.replace('JOBID',str(j))
            condor_file.write(line)
        condor_file.close()

        csh_file = open(dir[i]+"/"+prefix[i]+"_"+str(j)+".csh","w")
        for line in csh_templ_file:
            line=line.replace('CMSSWBASE',rel_base)
            line=line.replace('DIRECTORY',dir[i])
            line=line.replace('PREFIX',prefix[i])
            line=line.replace('JOBID',str(j))
            line=line.replace('FWLITEANALYZER',FWLiteAnalyzer)
            csh_file.write(line)
        csh_file.close()

        os.system('chmod u+x '+dir[i]+'/'+prefix[i]+'_'+str(j)+'.csh')
        os.system('cd '+dir[i]+'; condor_submit '+prefix[i]+'_'+str(j)+'.condor; cd -')
        j = j + 1
        nfiles = nfiles + files_per_job
        py_templ_file.close()
        condor_templ_file.close()
        csh_templ_file.close()

    #print  str(j-1)+' jobs submitted'
    
