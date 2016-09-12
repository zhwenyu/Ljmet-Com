import os, sys, getopt, datetime
execfile("/uscms_data/d3/varun/EOSSafeUtils.py")
Dir = sys.argv[1]

print; print 'Checking', Dir

try:
    opts, args = getopt.getopt(sys.argv[2:], "", ["verbose=", "resubmit=", "resub_num=", "doNoRoot="])
except getopt.GetoptError as err:
    print str(err)
    sys.exit(1)
    
verbose_level = 0
resubmit = '0'
resub_num = -2
doNoRoot = False
doNoLog = False

for o, a in opts:
	print o, a
	if o == '--verbose': verbose_level = int(a)
	if o == '--resubmit': resubmit = a
	if o == '--resub_num': resub_num = int(a)
	if o == '--doNoRoot': doNoRoot = bool(int(a))

folders = EOSlistSubdirs(Dir) #was folders = [x for x in os.walk(Dir).next()[1]] 

total_total = 0
total_succeeded = 0
total_error = 0
total_running = 0

no_log = 0
empty_log = 0
no_root = 0

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
time='%i_%i_%i'%(cTime.hour,cTime.minute,cTime.second)
tempdir = 'templogs_'+date+'_'+time
os.makedirs(tempdir)

for folder in folders:
 	#if folder != 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns': continue
    if 'Tprime' not in folder: continue
    if verbose_level > 0:  print; print folder
    if 'logfiles' not in folder: folder += '/logfiles/'
    files = [x for x in EOSlistdir(Dir+'/'+folder) if '.py' in x] #was os.listdir zzz
#Return a list containing the names of the entries in the directory given by path. 
#The list is in arbitrary order. It does not include the special entries '.' and '..' even if they are present in the directory.
	
    EOSlistdir(Dir+'/'+folder)  #was os.listdir
	
    resub_index = []
    count_total = 0
    for file in files:
        total_total+=1
        index = file[file.find('_')+1:file.find('.')]
        if '_' in index: index = index.split('_')[-1]
        count_total += 1

        if not EOSisfile(Dir+'/'+folder+'/'+file.replace('.py','.log')): #was os.path.isfile
            if verbose_level > 0: 
                print '\tNO LOG:',file,' and JobIndex:',index
            no_log += 1
            if resub_num == -1 or resub_num == 1: resub_index.append(index)
            continue

        current = EOSopen_via_temp(Dir + '/'+folder+'/'+file.replace('.py','.log'),'r') #open(...,'r')->EOSopen(...)
        good = False
        for line in current:
            if 'All cuts' in line: good = True
        current.close()
        if not good: 
            if verbose_level > 0: 
                print '\tEMPTY LOG:',file,' and JobIndex:',index
            empty_log+=1
            if resub_num == -1 or resub_num == 0:resub_index.append(index)
            continue

        if not EOSisfile(Dir+'/'+folder.replace('logfiles','')+file.replace('.py','.root')): #was os.path.isfile
            if verbose_level > 0: 
                print '\tNO ROOT:',file,' and JobIndex:',index
            no_root+=1
            if resub_num == -1 or resub_num == 2: resub_index.append(index)
            continue

    if resub_index != []: print 'RESUBS:', resub_index
    if resubmit != '1': continue

    indexind = 0
    rundir = os.getcwd()
    for index in resub_index:

 #       if 'ST' in folder and '7' in index: continue

        print '-----------------------------------------------'
        cleanpath = cleanEOSpath(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py')
        print 'resubmitting',cleanpath
        os.system('xrdcp -fs root://cmseos.fnal.gov/'+cleanpath+' '+tempdir+'/')

        tempfile = tempdir+'/'+folder.replace('/logfiles/','') + '_' + index + '.py'
        
        f = open(tempfile,'rU') #EOSopen_via_temp(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py', 'rU') #open(...,'r')->EOSopen(...)
        ConfigLines = f.readlines()
        f.close()
            #EOSopen_via_temp(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py','w') as fout:  #open(...,'r')->EOSopen(...) zzz
        with open(tempfile,'w') as fout:
            for line in ConfigLines:
                if line.startswith('relBase    = str('): 
                    print 'found the line to replace'
                    fout.write('relBase    = str(\'CONDOR_RELBASE\')\n')
                else: fout.write(line)
        os.system('xrdcp -f '+tempfile+' root://cmseos.fnal.gov/'+cleanpath)
        EOSrm(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.log') #os.system('rm '...) -> EOSrm(...)
        os.chdir(Dir + '/' + folder)
        print 'moved to',os.getcwd()
        jdlname = folder.replace('/logfiles/','') + '_' + index + '.jdl'
        print 'submitting',jdlname
        os.system('condor_submit ' + jdlname) #condor_submit is fixed to be eos safe.
        indexind+=1
        os.chdir(rundir)

print
print 'TOTAL JOBS: ', total_total
print 'NO LOG:', no_log
print 'EMPTY LOG:', empty_log
print 'NO ROOT:', no_root
print 'DONE:', total_total - no_log - empty_log - no_root
	
