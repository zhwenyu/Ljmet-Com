

import os
import sys

#print os.walk('/uscms_data/d3/rsyarif/Brown2018/TT_BB_SSDL_LJMet_2018data/LJMet102x_2lepTT_2018datasets_2018_11_29_rizki/').next()[1]

path = '/uscms_data/d3/rsyarif/Brown2018/TT_BB_SSDL_LJMet_2018data/LJMet102x_2lepTT_2018datasets_2018_11_29_rizki/'

try:
	samples = [x for x in os.walk(path).next()[1]]
except Exception as e:
	print e


oldMem = '8000'
newMem = '10000'

try:
	oldMem=sys.argv[2]
	newMem=sys.argv[3]
except:
	pass

total_changed=0
mem_fail_total=0

EditRequestMemory = 0 #False

try:
    EditRequestMemory = int(sys.argv[1])
except:
    pass

for sample in samples:
    #if 'TprimeTprime_M-1400' not in sample: continue
    #if 'DoubleMuon_RRF_v2' not in sample:continue

    print 'Checking:',sample

    files = [x for x in os.listdir(path+sample) if 'condor.log' in x]

    for f in files:

        with open(path+sample+'/'+f) as logf:
            line = logf.read()
            if 'exceeding requested memory' in line:

                print '\t',f
                print '\t\t--> contains "exceeding requested memory"'

                mem_fail_total+=1

                orig_f_lines = open(path+sample+'/'+f.replace('.log','')).readlines()
                for line in orig_f_lines:
                    if 'memory' in line:
                        print '\t\t\tCurrently in file: ',line

                if (EditRequestMemory==False): continue

                orig_f =  open(path+sample+'/'+f.replace('.log','')).read()
                if 'request_memory = '+oldMem in orig_f:
                    print '\t\t\t replacing request_memory',oldMem,'-->',newMem
                    orig_f = orig_f.replace('request_memory = '+oldMem,'request_memory = '+newMem)
                    print '\t\t\t writing to a NEW file!', path+sample+'/'+f.replace('.log','')
                    new_f = open(path+sample+'/'+f.replace('.log',''),'w') #THIS REWRITES EMPTY NEW FILE!!
                    new_f.write(orig_f)
                    new_f.close()

                    total_changed+=1


print '='*30
print total_changed,'files rewritten.'
print mem_fail_total,'jobs with memory fail.'
