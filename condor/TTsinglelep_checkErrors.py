import os, sys, getopt
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

folders = [x for x in os.walk(dir).next()[1]]

total_total = 0
total_succeeded = 0
total_error = 0
total_running = 0

no_log = 0
empty_log = 0

for folder in folders:
 	#if folder != 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns': continue
        #if 'Bprime' not in folder: continue
	if verbose_level > 0:  print; print folder

	files = [x for x in os.listdir(dir+'/'+folder) if '.jdl' in x]
	
	os.listdir(dir+'/'+folder)
	
	resub_index = []
	count_total = 0
	for file in files:
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
print 'DONE:', total_total - no_log - empty_log
