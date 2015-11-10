#!/usr/bin/python

import os
import re
import fileinput

files_per_job = 1

rel_base = os.environ['CMSSW_BASE']

outdir = '/eos/uscms/store/user/lpctlbsm/clint/Spring15/25ns/FakeRate/PreApproval/'

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
            #'TTJets',
            'TTbar-powheg',
            'TTbar-powheg_ext',
            #'QCD_Pt-20to30_EMEnriched',
            #'QCD_Pt-30to50_EMEnriched',
            #'QCD_Pt-50to80_EMEnriched',
            #'QCD_Pt-80to120_EMEnriched',
            #'QCD_Pt-120to170_EMEnriched',
            #'QCD_Pt-170to300_EMEnriched',
            #'QCD_Pt-300toInf_EMEnriched',
            #'QCD_Pt-20to30_MuEnrichedPt5',
            #'QCD_Pt-30to50_MuEnrichedPt5',
            #'QCD_Pt-50to80_MuEnrichedPt5',
            #'QCD_Pt-80to120_MuEnrichedPt5',
            #'QCD_Pt-120to170_MuEnrichedPt5',
            #'QCD_Pt-170to300_MuEnrichedPt5',
            #'QCD_Pt-300to470_MuEnrichedPt5',
            #'QCD_Pt-470to600_MuEnrichedPt5',
            #'QCD_Pt-600to800_MuEnrichedPt5',
            #'QCD_Pt-800to1000_MuEnrichedPt5',
            #'QCD_Pt-1000toInf_MuEnrichedPt5',
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
#    'Samples_Spring15MC/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74_Asympt25ns_MCRUN2_74_V9-v1.txt',
    'Samples_Spring15MC/TT_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
    'Samples_Spring15MC/TT_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3.txt',
#    'Samples_Spring15MC/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1.txt',
#    'Samples_Spring15MC/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
#    'Samples_Spring15MC/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2.txt',
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

        py_templ_file = open(rel_base+"/src/LJMet/Com/condor/Dilepton_Spring15MC_FakeRate_python.templ")
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
    
