from odbAccess import *
import sys

# script for extracting KKS actuator data

s=sys.argv[1]

odb = openOdb(path='%s'%(s))

time_list = []
ctf1_post_list = []
ctm1_post_list = []



fw=open("JOINT_LOADS.results",'w')
fw.write("TIME(sec) TF_CTF1 TF_CTM1   \n")

for elmset in odb.steps['Step-3'].historyRegions.keys():

############################### QUADRICEPS RESULTS ################################
   # quad loads
   if (elmset == 'Element PART-1-1.100000'):   # TF sensor
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data)
       for i in range(numData):
           time_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][0])
           ctf1_post_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
           ctm1_post_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTM1'].data[i][1])
           

for i in range(len(time_list)):

    fw.write(" %9.3f %9.3f %9.3f   \n"%(time_list[i],ctf1_post_list[i],ctm1_post_list[i]))

odb.close()
