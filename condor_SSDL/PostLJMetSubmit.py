#!/usr/bin/python


import os
import re
import fileinput
import commands

rel_base = os.environ['CMSSW_BASE']
cmssw = 'CMSSW_9_4_11'
date = 'Dec14'
folder  = 'LJMet94x_2lepTT_2017datasets_2018_11_18_rizki'

locdir = 'hadd_jobs_'+folder
basedir = '/store/group/lpcljm/'
indir = 'root://cmseos.fnal.gov/'+basedir+'/'+folder+'/'
outdir = 'root://cmseos.fnal.gov/'+basedir+'/'+folder+'_hadds/'
eosdir = '/store/group/lpcljm/'+folder+'_hadds/'

#################################################
### Names to give to your output root files
#################################################

samples = [
    'DoubleEG_RRB',
    'DoubleEG_RRC',
    'DoubleEG_RRD',
    'DoubleEG_RRE',
    'DoubleEG_RRF_v2',
    'MuonEG_RRB',
    'MuonEG_RRC',
    'MuonEG_RRD',
    'MuonEG_RRE',
    'MuonEG_RRF_v2',
    'DoubleMuon_RRB',
    'DoubleMuon_RRC',
    'DoubleMuon_RRD',
    'DoubleMuon_RRE',
    'DoubleMuon_RRF_v2',

   'TTW',
   'TTZ',
   'TTH',
   'TTTT',

   'WWW',
   'WWZ',
   'WZZ',
   'ZZZ',
   
   'WpWp',
   'WZ',
   'ZZ',
   
#    'TprimeTprime_M-1000',
   'TprimeTprime_M-1100',
   'TprimeTprime_M-1200',
   'TprimeTprime_M-1300',
   'TprimeTprime_M-1400',
   'TprimeTprime_M-1500',
   'TprimeTprime_M-1600',
   'TprimeTprime_M-1700',
   'TprimeTprime_M-1800',
   
]

#samples = []


datalist = [
]

systlist = [
# '_JESDOWN',
# '_JESUP',
# '_JERUP',
# '_JERDOWN',
'NOM'
]

#samples = ['Wprime3100Right','Wprime3200Right','Wprime3300Right','Wprime3400Right','Wprime3500Right','Wprime3600Right','Wprime3700Right','Wprime3800Right','Wprime3900Right','Wprime4000Right']

### Write the files you wish to run over for each job    

#make local directory
os.system('mkdir -p  %s' %locdir)
#make eos directory
os.system('eos root://cmseos.fnal.gov mkdir -p  %s' %eosdir)

for i in range(len(samples)):
    for sys in systlist:
        if(sys!='NOM'): continue
        condor_templ_file = open(rel_base+"/src/LJMet/Com/condor_SSDL/PostLJMetcondor.templ")
        csh_templ_file    = open(rel_base+"/src/LJMet/Com/condor_SSDL/PostLJMetcsh.templ")
    
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
				line=line.replace('HADD', '''hadd -f '''+samples[i]+'''_tmp.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+folder+'''/'''+samples[i]+'''/ | grep ".root"`''')
				line=line.replace('XRDCP', 'xrdcp -f '+samples[i]+'_tmp.root root://cmseos.fnal.gov/'+eosdir+'/'+samples[i]+'.root')
            else:
				line=line.replace('HADD', '''hadd -f '''+samples[i]+sys+'''_tmp.root `xrdfs root://cmseos.fnal.gov ls -u '''+basedir+'''/'''+folder+'''/'''+samples[i]+sys+'''/ | grep ".root"`''')
				line=line.replace('XRDCP', 'xrdcp -f '+samples[i]+sys+'_tmp.root root://cmseos.fnal.gov/'+eosdir+'/'+samples[i]+sys+'.root')
            csh_file.write(line)
        csh_file.close()
    
    
        if (sys=='NOM'): os.system('chmod u+x '+locdir+'/'+samples[i]+'.csh')
        else: os.system('chmod u+x '+locdir+'/'+samples[i]+sys+'.csh')
        print 'condor file is: '+localcondor
        if (os.path.exists('%s.log'  % localcondor)):
            os.system('rm %s.log' % localcondor)
        os.system('echo condor_submit %s' % localcondor)
        os.system('condor_submit %s' % localcondor)

    
        condor_templ_file.close()
        csh_templ_file.close()
    
