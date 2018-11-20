import os, sys, getopt
from eos_utils.EOSSafeUtils import *  # <--> file from "/uscms_data/d3/varun/EOSSafeUtils.py"

###Syntax:
#          CheckErrors.py /uscms_data/d3/path/to/logs/nominal/ --verbose 1 (prints) --resubmit 1 (if you want to resubmit) --resub_num -1 (resubmits all the fails it looks for)
#

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

rootdir = '/store/group/lpcljm/'+dir.split('/')[-2]+'/'
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
err_fail = 0
copy_fail = 0
mem_fail = 0 ; mem_fail_dict = {}
kill_fail = 0
All_cuts = 0

root_not_found = {}

for folder in folders:
# 	if 'TprimeTprime_M-1400' not in folder: continue
# 	if 'DoubleMuon_RRB' not in folder: continue
# 	if 'ZZ' not in folder: continue
	if verbose_level > 0:  
		print
		print folder
        rootfiles = EOSlist_root_files(rootdir+folder)
        total_roots += len(rootfiles)
        
        mem_fail_dict[folder] = [] #initialize list, for saving memory failed jobs
        root_not_found[folder] = [] #initialize list, for finished jobs with no root files


	files = [x for x in os.listdir(dir+'/'+folder) if '.condor.log' in x]
	
	os.listdir(dir+'/'+folder)
	
	resub_index = []
	count_total = 0
	for file in files:
		#print file
		total_total+=1
		index = file[file.find('_')+1:file.find('.')]
		if '_' in index: index = index.split('_')[-1]
		count_total += 1
	
		#Check if stdout exists
		if (os.path.exists(dir+'/'+folder+'/'+file.replace('.condor.log','.stdout'))==False):
	
			print '\tNO STDOUT:',file.replace('.condor.log','.stdout'),' and JobIndex:',index
			no_log += 1
			#if resub_num == -1 or resub_num == 1: resub_index.append(index)

			#Check why does it not have stdout?
			try:
				current = open(dir + '/'+folder+'/'+file.replace('.condor.log','.condor.log'),'r')
				overmem = False
				killed = False
				term = False
				iline = 0
				overmem_index = 0
				kill_index = 0
				term_index = 0
				for line in current:
					if 'SYSTEM_PERIODIC_REMOVE' in line: 
										overmem = True
										overmem_index = iline
					elif 'condor_rm' in line: 
						killed = True
						kill_index = iline
					elif 'Normal termination (return value 0)' in line: 
						term = True
						term_index = iline
					iline += 1
				if overmem and overmem_index > term_index: 
					if verbose_level > 0: 
						print '\t\tMEM FAIL:',file.replace('.condor.log','.condor.log'),' and JobIndex:',index
					mem_fail+=1
					if resub_num == -1 or resub_num == 3:resub_index.append(index)
					mem_fail_dict[folder].append(index)
					continue
				if killed and kill_index > term_index: 
					if verbose_level > 0: 
						print '\t\tKILL_FAIL:',file.replace('.condor.log','.condor.log'),' and JobIndex:',index
					kill_fail+=1
					if resub_num == -1 or resub_num == 4:resub_index.append(index)
					continue
			except:
				pass
			

		else:
			current = open(dir + '/'+folder+'/'+file.replace('.condor.log','.stdout'),'r')
			good = False
			for line in current:
				if 'All cuts' in line: 
					good = True
					All_cuts +=1
					
					#check if root file exists:
					if file.replace('.condor.log','.root') not in rootfiles:
						print "\tROOT file does not exist for JobIndex:", index,'!!'
						root_not_found[folder].append(index)

						#check why doesn't root file exist?
						stderr_ = open(dir + '/'+folder+'/'+file.replace('.condor.log','.stderr'),'r')
						for l in stderr_:
							if 'Error' in l:
								if verbose_level > 0: 
									print '\t\t"Error" in STDERR:',file.replace('.condor.log','.stderr'),' and JobIndex:',index
								err_fail+=1
						stderr_.close()

						if resub_num == -1 or resub_num == 0:resub_index.append(index)
						

			#Check why does it not have "All Cuts" in stdout?
			if not good: 
				if verbose_level > 0: 
					print '\tNO "All Cuts" in STDOUT:',file.replace('.condor.log','.stdout'),' and JobIndex:',index
					empty_log+=1

					try:
						current = open(dir + '/'+folder+'/'+file.replace('.condor.log','.stderr'),'r')
						good = True
						for line in current:
							if 'Error' in line: good = False
						if not good: 
							if verbose_level > 0: 
								print '\t\t"Error" in STDERR:',file.replace('.condor.log','.stderr'),' and JobIndex:',index
							err_fail+=1
					except:
						pass

					if resub_num == -1 or resub_num == 0:resub_index.append(index)

					continue
# 		except:
# 			pass
# 
# 		try:
			current = open(dir + '/'+folder+'/'+file.replace('.condor.log','.stderr'),'r')
			good = True
			for line in current:
				if 'failure' in line: good = False
			if not good: 
				if verbose_level > 0: 
					print '\t"failure" in STDERR:',file.replace('.condor.log','.stderr'),' and JobIndex:',index
				copy_fail+=1
				if resub_num == -1 or resub_num == 2:resub_index.append(index)
				continue
# 		except:
# 			pass

		
	if resub_index != []: print 'RESUBS:', resub_index
	if resubmit != '1': continue

	indexind = 0
	for index in resub_index:
		os.chdir(dir + '/' + folder)

		os.system('rm -v ' + dir + '/' + folder + '/' + folder + '_' + index + '.condor.log')
		try:
			os.system('rm -v ' + dir + '/' + folder + '/' + folder + '_' + index + '.stdout')
			os.system('rm -v ' + dir + '/' + folder + '/' + folder + '_' + index + '.stderr')
		except:
			pass
		os.system('condor_submit ' + dir + '/' + folder + '/' + folder + '_' + index + '.condor')

		indexind+=1
	
print	
print '=' * 50
print
print 'TOTAL JOBS: ', total_total
print
print 'ROOT files found:', total_roots
print '\tSTDOUT prints "All Cuts":', All_cuts
print '\tSTDOUT prints no "All Cuts":', empty_log
print '\t\t"Error" in STDERR:', err_fail
#print 'DONE:', total_total - no_log - empty_log - copy_fail - mem_fail - kill_fail
print
print 'NO STDOUT:', no_log
print '\tMEMORY FAIL:', mem_fail
# for sample in mem_fail_dict:
# 	if mem_fail_dict[sample]==[]: continue
# 	print '\t\t ',sample, mem_fail_dict[sample]
print '\tKILL_FAIL:', kill_fail
print '\tProbably still running:', no_log - mem_fail - kill_fail
print
print '"failure" in STDERR:', copy_fail
print
print 'DONE (TOTAL JOBS - NO STDOUT - NO "All Cuts" - "failure" FAILS):', total_total - no_log - empty_log - copy_fail
print
print '( "All Cuts" + NO "All Cuts" + NO STDOUT + "failure" FAILS):', All_cuts + no_log + empty_log + copy_fail
print

