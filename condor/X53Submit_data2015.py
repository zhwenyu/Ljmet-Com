#!/usr/bin/python

import os
import re
import fileinput

files_per_job = 1

rel_base = os.environ['CMSSW_BASE']

###########################################
### Where to save your output root files to
###########################################
outdir = '/uscms_data/d3/drankin/13TeV/Samples/13Aug_NoIso/'

### What is the name of your FWLite Analyzer
FWLiteAnalyzer = 'ljmet'

### Selection to run
SELECTOR = "'singleLepSelector'"

### Is this sample MC?
ISMCSAMPLE = 'False'

###Is WJets?
ISWJETS = 'False'

###53x JEC?
DONEWJEC = 'False'

###Is TB?
ISTB = 'False'
ISTT = 'False'

### Use best top for neutrino pz solution?
USEBESTTOP = 'True'

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

prefix = [
'SingleMuon_Run2015B_PromptReco',
'SingleMuon_Run2015B_17Jul2015',
'SingleElectron_Run2015B_PromptReco',
'SingleElectron_Run2015B_17Jul2015',
]

dir = [
outdir+'SingleMuon_Run2015B_PromptReco',
outdir+'SingleMuon_Run2015B_17Jul2015',
outdir+'SingleElectron_Run2015B_PromptReco',
outdir+'SingleElectron_Run2015B_17Jul2015',
]


################################################
### Where is your list of root files to run over
################################################
list = [ 
'Samples_2015/SingleMuon_Run2015B_PromptReco.txt',
'Samples_2015/SingleMuon_Run2015B_17Jul2015.txt',
'Samples_2015/SingleElectron_Run2015B_PromptReco.txt',
'Samples_2015/SingleElectron_Run2015B_17Jul2015.txt',
] 

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
                f_name=re.search('.+\'(.+\.root)',line)
                #result=result+'                 \'' + f_name.group(1)+'\',\n'
                #result=result+'                 \'dcap:///pnfs/cms/WAX/11' + f_name.group(1)+'\',\n'
                result=result+'                 \'root://cmsxrootd.fnal.gov/' + f_name.group(1)+'\',\n'
    file_list.close()
    #result = result + '                 )\n'
    return result


print str(files_per_job)+' files per job...'

for i in range(len(prefix)):

    j = 1
    nfiles = 1

    #MYJSON = "''"
    #MYJSON = "'../data/json/json_DCSONLY_Run2015B.txt'"
    MYJSON = "'../data/json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'"
    #if (prefix[i].find('13Jul2012')>0):
    #    MYJSON = "'../data/json/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt'"
    #if (prefix[i].find('06Aug2012')>0):
    #    MYJSON = "'../data/json/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt'"
    #if (prefix[i].find('24Aug2012')>0):
    #    MYJSON = "'../data/json/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'"
    #if (prefix[i].find('Prompt2012C')>0):
    #    MYJSON = "'../data/json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_198941-203709.txt'"
    #if (prefix[i].find('Prompt2012D')>0):
    #    MYJSON = "'../data/json/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON_203894-208686.txt'"
    #if (prefix[i].find('11Dec2012')>0):
    #    MYJSON = "'../data/json/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt'"

    FLAGTAG = 'TriggerResults::RECO'
    if (prefix[i].find('17Jul2015')>0):
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

        #if (prefix[i].startswith('SingleElectron') and j <= 24):
        #    j = j + 1
        #    nfiles = nfiles + files_per_job
        #    continue
        #if (prefix[i].startswith('SingleMuon') and j <= 16):
        #    j = j + 1
        #    nfiles = nfiles + files_per_job
        #    continue

        #py_templ_file = open(rel_base+"/src/LJMet/Com/condor/Wprimepython.templ")
        py_templ_file = open(rel_base+"/src/LJMet/Com/condor/Wprimepython.templ.data")
        condor_templ_file = open(rel_base+"/src/LJMet/Com/condor/Wprimecondor.templ")
        csh_templ_file    = open(rel_base+"/src/LJMet/Com/condor/Wprimecsh.templ")

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

    print  str(j-1)+' jobs submitted'
    
