import os, sys, getopt
execfile("/uscms_data/d2/abarker/area2/CMSSW_7_4_14/src/LJMet/Com/condor_2lep_120715/EOSSafeUtils.py")
Dir = sys.argv[1]

print; print 'Checking', Dir

try:
    opts, args = getopt.getopt(sys.argv[2:], "", ["verbose=", "resubmit="])
except getopt.GetoptError as err:
    print str(err)
    sys.exit(1)
    
verbose_level = 0
resubmit = '0'
doNoLog = False

for o, a in opts:
	print o, a
	if o == '--verbose': verbose_level = int(a)
	if o == '--resubmit': resubmit = a

folders = EOSlistSubdirs(Dir) #was folders = [x for x in os.walk(Dir).next()[1]] 

total_total = 0
total_succeeded = 0
total_error = 0
total_running = 0

no_log = 0
for folder in folders:
#    if 'Tprime' not in folder: continue
#    if '1000' not in folder: continue
    if verbose_level > 0:  print; print folder
    if 'logfiles' not in folder: folder += '/logfiles'
    files = [x for x in EOSlistdir(Dir+'/'+folder) if 'jdl' in x] #was os.listdir zzz
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
        
	try:
            namestring = Dir+folder+'/'+file.replace('.jdl','.condor')
            nameshort = namestring[10:]
            os.system('xrdcp -s root://cmseos.fnal.gov/'+nameshort+' ./.eostemp')
            current = open('./.eostemp','r')
            good = False
            for line in current:
                if 'Job submitted' in line: good = True
            current.close()
            os.system('rm ./.eostemp')
            if not good: 
                os.system('xrdcp -s root://cmseos.fnal.gov/'+nameshort.replace('.condor','.py')+' ./.eostemp')
                current = open('./.eostemp','r')
                relbaseok = False
                for line in current:
                    if 'CONDOR_RELBASE' in line: relbaseok = True
                current.close()
                os.system('rm ./.eostemp')

                if verbose_level > 0: 
                    print '\tNot submitted:',file,' and JobIndex:',index,' relbase ok?',relbaseok
	       	no_log+=1
                resub_index.append(index)
                continue
	except:
            print 'that check failed'
            pass

    if resub_index != []: print 'RESUBS:', resub_index
    if resubmit != '1': continue
        
    indexind = 0
    rundir = os.getcwd()
    for index in resub_index:
        namestring = folder.replace('/logfiles','') + '_' + index + '.jdl'
        dirstring = Dir + folder + '/'
        os.chdir(dirstring)
        os.system('condor_submit ' + namestring)
        indexind+=1
        os.chdir(rundir)
	
	

print
print 'TOTAL JOBS: ', total_total
print 'NO LOG:', no_log
print 'DONE:', total_total - no_log
	
