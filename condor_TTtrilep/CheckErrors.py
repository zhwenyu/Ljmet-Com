import os, sys, getopt
execfile("/uscms_data/d3/varun/EOSSafeUtils.py")

## CheckErrors.py /uscms_data/d3/path/to/logs/nominal/ --verbose 1 (prints) --resubmit 1 (if you want to resubmit) --resub_num -1 (resubmits all the fails it looks for)

dir = sys.argv[1]

print; print 'Checking', dir

try:
    opts, args = getopt.getopt(sys.argv[2:], "", ["verbose=", "resubmit=", "resub_num="])
except getopt.GetoptError as err:
    print str(err)
    sys.exit(1)
    
verbose_level = 0
resubmit = '0'
resub_num = -2
doNoLog = False

for o, a in opts:
	print o, a
	if o == '--verbose': verbose_level = int(a)
	if o == '--resubmit': resubmit = a
	if o == '--resub_num': resub_num = int(a)

rootdir = '/eos/uscms/store/user/lpcljm/'+dir.split('/')[-3]+'/'+dir.split('/')[-2]+'/'
rootdir = rootdir.replace('_logs','')
print 'checking ROOT files in:',rootdir
folders = [x for x in os.walk(dir).next()[1]]

total_total = 0
total_succeeded = 0
total_error = 0
total_running = 0
total_roots = 0

no_log = 0
empty_log = 0
copy_fail = 0

for folder in folders:
        #if 'DoubleEG' not in folder: continue 
        #if folder != 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns': continue
        #if 'Charged' not in folder and 'Single' not in folder: continue
	if verbose_level > 0:  print; print folder

        rootfiles = EOSlist_root_files(rootdir+folder)
        total_roots += len(rootfiles)

	files = [x for x in os.listdir(dir+'/'+folder) if '.jdl' in x]
	
	os.listdir(dir+'/'+folder)
	
	resub_index = []
	count_total = 0
	for file in files:
		#print file
		total_total+=1
		index = file[file.find('_')+1:file.find('.')]
		if '_' in index: index = index.split('_')[-1]
		count_total += 1
	
		try:
			current = open(dir + '/'+folder+'/'+file.replace('.jdl','.log'),'r')
			good = False
			for line in current:
				if 'All cuts' in line: good = True
			if not good: 
				if verbose_level > 0: 
					print '\tEMPTY LOG:',file,' and JobIndex:',index
				empty_log+=1
				if resub_num == -1 or resub_num == 0:resub_index.append(index)
				continue
		except:
			pass

		try:
			current = open(dir + '/'+folder+'/'+file.replace('.jdl','.out'),'r')
			good = True
			for line in current:
				if 'failure' in line: good = False
			if not good: 
				if verbose_level > 0: 
					print '\tXRDCP FAIL:',file,' and JobIndex:',index
				copy_fail+=1
				if resub_num == -1 or resub_num == 2:resub_index.append(index)
				continue
		except:
			pass
		
		try:
			if not os.path.isfile(dir+'/'+folder+'/'+file.replace('.jdl','.log')): 
				if verbose_level > 0: 
					print '\tNO LOG:',file,' and JobIndex:',index
				no_log += 1
				if resub_num == -1 or resub_num == 1: resub_index.append(index)
				continue
		except: pass

	if resub_index != []: print 'RESUBS:', resub_index
	if resubmit != '1': continue

	indexind = 0
	for index in resub_index:
		os.chdir(dir + '/' + folder)

                f = open(dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py', 'rU')
                ConfigLines = f.readlines()
                f.close()
                with open(dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py','w') as fout:
                    for line in ConfigLines:
			if line.startswith('relBase    = str('): fout.write('relBase    = str(\'CONDOR_RELBASE\')\n')
			else: fout.write(line)
		os.system('rm ' + dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.log')
		os.system('condor_submit ' + dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.jdl')
		indexind+=1
	
	

print
print 'TOTAL JOBS: ', total_total
print 'NO LOG:', no_log
print 'EMPTY LOG:', empty_log
print 'XRDCP FAIL:', copy_fail
print 'ROOT files:', total_roots
print 'DONE:', total_total - no_log - empty_log - copy_fail
