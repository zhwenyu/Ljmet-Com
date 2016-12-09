#!/usr/bin/python


import os
import re
import fileinput
import commands

rel_base = os.environ['CMSSW_BASE']
cmssw = 'CMSSW_8_0_22'
date = 'Nov16'
locdir = date
basedir = '/store/user/clint/Run2016'
indir = 'root://cmseos.fnal.gov/'+basedir+'/'+date+'/'
outdir = 'root://cmseos.fnal.gov/'+basedir+'/'+date+'/'
eosdir = ''+basedir+'/'+date+'/'

#################################################
### Names to give to your output root files
#################################################

samples = [
#    'Data',
#    'DoubleEG-Run2016B-23Sep2016-v3',
#    'DoubleEG-Run2016C-23Sep2016-v1',
#    'DoubleEG-Run2016D-23Sep2016-v1',
#    'DoubleEG-Run2016E-23Sep2016-v1',
#    'DoubleEG-Run2016F-23Sep2016-v1',
#    'DoubleEG-Run2016G-23Sep2016-v1',
#    'DoubleEG-Run2016H-PromptReco-v2',
#    'DoubleEG-Run2016H-PromptReco-v3',
#    'DoubleMuon-Run2016B-23Sep2016-v3',
#    'DoubleMuon-Run2016C-23Sep2016-v1',
    'DoubleMuon-Run2016D-23Sep2016-v1',
#    'DoubleMuon-Run2016E-23Sep2016-v1',
#    'DoubleMuon-Run2016F-23Sep2016-v1',
#    'DoubleMuon-Run2016G-23Sep2016-v1',
#    'DoubleMuon-Run2016H-PromptReco-v2',
#    'DoubleMuon-Run2016H-PromptReco-v3',
#    'MuonEG-Run2016B-23Sep2016-v3',
#    'MuonEG-Run2016C-23Sep2016-v1',
#    'MuonEG-Run2016D-23Sep2016-v1',
#    'MuonEG-Run2016E-23Sep2016-v1',
#    'MuonEG-Run2016F-23Sep2016-v1',
#    'MuonEG-Run2016G-23Sep2016-v1',
#    'MuonEG-Run2016H-PromptReco-v2',
#    'MuonEG-Run2016H-PromptReco-v3',
#    'TTW',
#    'TTZ',
#    'TTH',
]

rmasses = [
'1000',
'1100',
'1200',
'1300',
'1400',
'1500',
'1600',
'1700',
'1800',
'1900',
'2000',
'2100',
'2200',
'2300',
'2400',
'2500',
'2600',
'2700',
'2800',
'2900',
'3000',
]

#for i in rmasses:
#    samples.extend(['Wprime'+i+'Right'])

datalist = [
]

systlist = [
#'_JESDOWN',
#'_JESUP',
#'_JERUP',
#'_JERDOWN',
'NOM'
]

#samples = ['Data']

### Write the files you wish to run over for each job    

#make local directory
#os.system('mkdir -p  %s' %locdir)
#make eos directory
#os.system('eos root://cmseos.fnal.gov mkdir -p  %s' %eosdir)

for i in range(len(samples)):
    for sys in systlist:
#        if (sys!='NOM' and (samples[i].startswith('QCD') or samples[i].startswith('TTbar_scale') or samples[i]=='Data')): continue
        if(sys!='NOM'): continue
        condor_templ_file = open(rel_base+"/src/LJMet/Com/condor/PostLJMetcondor.templ")
        csh_templ_file    = open(rel_base+"/src/LJMet/Com/condor/PostLJMetcsh.templ")
    
        localcondor = locdir+'/'+samples[i]+sys+".condor"
        if (sys=='NOM'): localcondor = locdir+'/'+samples[i]+".condor"
        condor_file = open(localcondor,"w")
        for line in condor_templ_file:
            line=line.replace('DIRECTORY',locdir)
            if (sys=='NOM'): line=line.replace('PREFIX',samples[i])
            else: line=line.replace('PREFIX',samples[i]+sys)
            condor_file.write(line)
        condor_file.close()
    
        localcsh=locdir+'/'+samples[i]+sys+".csh"
        if (sys=='NOM'): localcsh=locdir+'/'+samples[i]+".csh"
        csh_file = open(localcsh,"w")
        for line in csh_templ_file:
            line=line.replace('CMSSWVERSION',cmssw)
            if (sys=='NOM'):
                if (samples[i] == 'ZJets_M50'):
                    line=line.replace('HADD', '''hadd -f ZJets_tmp.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+date+'''/ZJets_M50/ | grep ".root"`''')
                    line=line.replace('XRDCP', 'xrdcp -f ZJets_tmp.root root://cmseos.fnal.gov/'+eosdir+'/ZJets.root')
                elif (samples[i].find('2016')!=-1):
                    ida = 0
                    tmpstr = ''                    
                    gooddata = 1
                    firstrun = commands.getstatusoutput('eos root://cmseos.fnal.gov ls '+basedir+'/'+date+'/'+samples[i]+'/ | grep "_0.root"')
                    if (firstrun[1]==''):
                        numfile = commands.getstatusoutput('ls '+date+'/'+samples[i]+'/'+samples[i]+'*.stdout | wc -l')
                        for j in range(1,int(numfile[1])+1):
                            final = commands.getstatusoutput('grep -c \"All cuts          0\" '+date+'/'+samples[i]+'/'+samples[i]+'_'+str(j)+'.stdout')
                            if (int(final[1])!=1):
                                gooddata = j
                                break
                        tmpstr = tmpstr+'xrdcp -f root://cmseos.fnal.gov/'+basedir+'/'+date+'/'+samples[i]+'/'+samples[i]+'_'+str(gooddata)+'.root root://cmseos.fnal.gov/'+basedir+'/'+date+'/'+samples[i]+'/'+samples[i]+'_0.root\n'
                        tmpstr = tmpstr+'eosrm '+basedir+'/'+date+'/'+samples[i]+'/'+samples[i]+'_'+str(gooddata)+'.root\n'
                    tmpstr = tmpstr+'''hadd -f tmp_'''+str(ida)+'''.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+date+'''/'''+samples[i]+'''/ | grep ".root"`\n'''
                    ida+=1
                    line=line.replace('HADD', tmpstr+'hadd -f tmp.root tmp_*.root')
                    line=line.replace('XRDCP', 'xrdcp -f tmp.root root://cmseos.fnal.gov/'+eosdir+'ljmet_trees/ljmet_Data_'+samples[i]+'.root')
                else:
                    line=line.replace('HADD', '''hadd -f '''+samples[i]+'''_tmp.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+date+'''/'''+samples[i]+'''/ | grep ".root"`''')
                    line=line.replace('XRDCP', 'xrdcp -f '+samples[i]+'_tmp.root root://cmseos.fnal.gov/'+eosdir+'/'+samples[i]+'.root')
            else:
                if (samples[i] == 'ZJets_M50'):
                    line=line.replace('HADD', '''hadd -f ZJets'''+sys+'''_tmp.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+date+'''/ZJets_M50'''+sys+'''/ | grep ".root"`''')
                    line=line.replace('XRDCP', 'xrdcp -f ZJets'+sys+'_tmp.root root://cmseos.fnal.gov/'+eosdir+'/ZJets'+sys+'.root')
                else:
                    line=line.replace('HADD', '''hadd -f '''+samples[i]+sys+'''_tmp.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+date+'''/'''+samples[i]+sys+'''/ | grep ".root"`''')
                    line=line.replace('XRDCP', 'xrdcp -f '+samples[i]+sys+'_tmp.root root://cmseos.fnal.gov/'+eosdir+'/'+samples[i]+sys+'.root')
            csh_file.write(line)
        csh_file.close()
    
    
        if (sys=='NOM'): os.system('chmod u+x '+locdir+'/'+samples[i]+'.csh')
        #else: os.system('chmod u+x '+locdir+'/'+samples[i]+sys+'.csh')
        print 'condor file is: '+localcondor
        if (os.path.exists('%s.log'  % localcondor)):
            os.system('rm %s.log' % localcondor)
        os.system('condor_submit %s' % localcondor)
    
        condor_templ_file.close()
        csh_templ_file.close()
    
