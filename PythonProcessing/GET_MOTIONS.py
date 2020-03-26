from odbAccess import *
import sys

# script for extracting KKS actuator data

s=sys.argv[1]

odb = openOdb(path='%s'%(s))

time_list = []
tfap_list = []
tfie_list = []


fw=open("MOTIONS.results",'w')
fw.write("TIME(sec) TF_AP_TRANS TF_IE_ROT  \n")

for elmset in odb.steps['Step-3'].historyRegions.keys():

############################### QUADRICEPS RESULTS ################################
   # quad loads
   if (elmset == 'Element PART-1-1.151'):   # TF AP TRANSLATION
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       for i in range(numData):
           time_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][0])
           tfap_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])

   if (elmset == 'Element PART-1-1.150'):   # TF IE ROTATION
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CUR1'].data)
       for i in range(numData):
           tfie_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CUR1'].data[i][1])


for i in range(len(time_list)):

    fw.write(" %9.3f %9.3f %9.3f  \n"%(time_list[i],tfap_list[i],tfie_list[i]))

odb.close()
