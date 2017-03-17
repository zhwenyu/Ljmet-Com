import os, sys, getopt
import string
#This file has some general EOS safe utilities for eos file operations
# EOSpathExists(path)
# EOSisfile(path)
# EOSrm(path):
# EOSrmdir(path): 

xrd = 'eos root://cmseos.fnal.gov/'
#xrd = 'xrdfs root://cmseos.fnal.gov/' #doesn't work. 

def cleanEOSpath(path): 
	#if path starts with /eos/uscms remove it. 
	if string.find(path,'/eos/uscms',0,10) == 0:
		return path[10:]
	return path

def EOSpathExists(path): 
	#returns a bool true iff the path exists
	path = cleanEOSpath(path)
	return len(os.popen(xrd+' ls -d '+path).readlines()) == 1

def EOSisfile(path): 
	#returns a bool saying if this is a file. 
	path = cleanEOSpath(path)
	if not EOSpathExists(path):
		return False
	lines = os.popen(xrd+' ls -l '+path).readlines()
	if len(lines) != 1:
		return False
	else:
		return string.find(lines[0],'-',0,1) == 0 #entry should start with -

def EOSrm(path):
	#rm path, but EOS safe.
	path = cleanEOSpath(path)
	if EOSpathExists(path):
		os.system(xrd +" rm "+path)
	
def EOSrmdir(path): 
	#rmdir path, but EOS safe.
	path = cleanEOSpath(path)
	if EOSpathExists(path):
		os.system(xrd +" rmdir "+path)


def striplist(alist): 
	#takes a list of strings, returns a version of the list with 
	#whitespace stripped from all entries.
	ret = []
	for item in alist:
		ret.append(item.strip())
	return ret
	
def EOSlistdir(Dir): 
	#ls Dir, only eos safe.
	#does not check that this is a directory
	#returns a list of the contents of ls 
	Dir = cleanEOSpath(Dir)
	items = os.popen(xrd+' ls '+Dir).readlines() #they have a \n at the end 
	return striplist(items)

tempfile = 'eossafetemp'
def killtemp(): 
	#deletes the the temp file 
	if os.path.exists(tempfile):
		os.system("rm "+tempfile)

def EOSopen_via_temp(path): 
	#copies the file to the temp file, then opens that. 
	#be careful to not open more than one file at a time 
	#or else the first file will get deleted. 
	path = cleanEOSpath(path)
	killtemp()
	os.system("xrdcp root://cmseos.fnal.gov/"+path+" ./"+tempfile)
	return open(tempfile,'r')

def copytotemp(path): 
	#copies a file from eos to the local temp file
	killtemp()
	path = cleanEOSpath(path)
	os.system("xrdcp root://cmseos.fnal.gov/"+path+" ./"+tempfile)
	
	#I need to replace thei line:
	#folders = [x for x in os.walk(Dir).next()[1] ] 
	#this lists all the folders in Dir. 

def EOSlistSubdirs(Dir):
	#returns a list of all the subdirectories of Dir
	#built to substitutes for the line
	#folders = [x for x in os.walk(Dir).next()[1] ] 
	Dir = cleanEOSpath(Dir)
	ret = []
	if not EOSpathExists(Dir): 
		return ret
	for line in os.popen(xrd+' ls -l '+Dir).readlines():
		words = line.split()
		if len(words) < 9: #if not a full ls -l entry
			continue
		elif words[0][0]=='d': #if its a directory
			ret.append(words[8]) #append the directory name
	return ret

#///////////////////////////////////////////////////////
#///////////////////////////////////////////////////////
#///////////////////////////////////////////////////////
#///////////////////////////////////////////////////////
#///////////////////////////////////////////////////////

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
for folder in folders:
 	#if folder != 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns': continue
#        if 'TT_' not in folder: continue
	if 'Single' in folder: continue
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
	
		try:
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
		except:
			pass
		
		try:
			if not EOSisfile(Dir+'/'+folder+'/'+file.replace('.py','.condor')): #was os.path.isfile
				if verbose_level > 0: 
					print '\tNO LOG:',file,' and JobIndex:',index
				no_log += 1
				if resub_num == -1 or resub_num == 1: resub_index.append(index)
				continue
		except: pass
		if not EOSisfile(Dir+'/'+folder.replace('logfiles','')+file.replace('.py','.root')): #was os.path.isfile
			if verbose_level > 0: 
				print '\tNO ROOT:',file,' and JobIndex:',index
			no_root+=1
			if resub_num == -1 or resub_num == 2: resub_index.append(index)
			continue
	if resub_index != []: print 'RESUBS:', resub_index
	if resubmit != '1': continue

	indexind = 0
	for index in resub_index:
		#os.chdir(Dir + '/' + folder) #not needed?? zzz
		#print os.getcwd()
                f = EOSopen_via_temp(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py', 'rU') #open(...,'r')->EOSopen(...)
                ConfigLines = f.readlines()
                f.close()
                #with open(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py','w') as fout:
                with EOSopen_via_temp(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.py','w') as fout:  #open(...,'r')->EOSopen(...) zzz
                    for line in ConfigLines:
			if line.startswith('relBase    = str('): fout.write('relBase    = str(\'CONDOR_RELBASE\')\n')
#                        elif line.startswith('    \'LjetsTopoCalc\','): fout.write('    \'LjetsTopoCalc\',\'LjetsTopoCalcNew\',\n')
			else: fout.write(line)
		EOSrm(Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.log') #os.system('rm '...) -> EOSrm(...)
		os.system('condor_submit ' + Dir + '/' + folder + '/' + folder.replace('/logfiles/','') + '_' + index + '.jdl') #condor_submit is fixed to be eos safe.
		indexind+=1
	
	

print
print 'TOTAL JOBS: ', total_total
print 'NO LOG:', no_log
print 'EMPTY LOG:', empty_log
print 'NO ROOT:', no_root
print 'DONE:', total_total - no_log - empty_log - no_root
	
'''
try:
	current = EOSopen_via_temp(Dir + '/'+folder+'/'+file.replace('.py','.log'),'r')
except:
	count_running += 1
#            if verbose: print '\tRUNNING:',file
	if doNoRoot: resub_index.append(index)
	continue

good = False
for line in current:
	if 'All cuts' in line: good = True
	#os.path.isfile -> EOSisfile
#         print doNoRoot and (not EOSisfile(Dir + '/'+folder+'/'+file.replace('.py','.log')))
current.close()
if not good or (doNoRoot and not EOSisfile(Dir + '/'+folder+'/'+file.replace('.py','.log'))):
	count_failed += 1
	if verbose==3: print '\tFAILED:',file
	resub_index.append(index)

#     print resub_index
'''
'''
if resubmit != '1': continue

for index in resub_index:
	read_file = EOSopen(Dir + '/'+folder+'/' + folder + '.job','r') #open(...,'r')->EOSopen(...)
	writ_file = open(Dir + '/'+folder+'/' + folder + '_' + index + '.job', 'w')

	for line in read_file:
		line = line.replace('$(process)', index)
		if 'Queue' in line: line = 'Queue 1\n'

		writ_file.write(line)

	read_file.close()
	writ_file.close()

	#os.chdir(Dir + '/'+ folder +'/') #guessing these are unnecessary
	os.system('rm ' + folder + '_' + index + '.log')
	os.system('condor_submit ' + folder + '_' + index + '.job')
	#os.chdir('../../../') #guessing this is unnecessary
'''

# print
# print 'Total Jobs:', total_total
# print 'Total Running:', total_running
# print 'Total Succeeded:', total_succeeded
# print 'Total Failed:', total_error
