from odbAccess import *
import sys

# script for extracting KKS actuator data

s=sys.argv[1]

odb = openOdb(path='%s'%(s))

time_list = []
als_CU1_list = []
als_CTF1_list = []
pfl_CU1_list = []
pfl_CTF1_list = []
pcapl_CU1_list = []
pcapl_CTF1_list = []
pcapm_CU1_list = []
pcapm_CTF1_list = []

lcla_CU1_list = []
lcla_CTF1_list = []
lclm_CU1_list = []
lclm_CTF1_list = []
lclp_CU1_list = []
lclp_CTF1_list = []

smcla_CU1_list = []
smcla_CTF1_list = []
smclm_CU1_list = []
smclm_CTF1_list = []
smclp_CU1_list = []
smclp_CTF1_list = []




fw=open("FORCE_LIGS.results",'w')
fw.write("TIME(sec) als_CU1 als_CTF1 pfl_CU1 pfl_CTF1 pcapl_CU1 pcapl_CTF1 pcapm_CU1 pcapm_CTF1 lcla_CU1 lcla_CTF1 lclm_CU1 lclm_CTF1 lclp_CU1 lclp_CTF1 smcla_CU1 smcla_CTF1 smclm_CU1 smclm_CTF1 smclp_CU1 smclp_CTF1 \n")

for elmset in odb.steps['Step-3'].historyRegions.keys():

############################### QUADRICEPS RESULTS ################################
   # quad loads
   if (elmset == 'Element PART-1-1.18000'):   # ALS TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           time_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][0])
           als_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           als_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])


   if (elmset == 'Element PART-1-1.43000'):   # PFL TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           pfl_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           pfl_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])


   if (elmset == 'Element PART-1-1.23000'):   # PCAPL TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           pcapl_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           pcapl_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
           
           
   if (elmset == 'Element PART-1-1.28000'):   # PCAPM TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           pcapm_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           pcapm_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])


   if (elmset == 'Element PART-1-1.10422'):   # LCLA TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           lcla_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           lcla_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
           
   if (elmset == 'Element PART-1-1.10424'):   # LCLM TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           lclm_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           lclm_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
          
   if (elmset == 'Element PART-1-1.10426'):   # LCLP TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           lclp_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           lclp_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])



   if (elmset == 'Element PART-1-1.33002'):   # SMCLA TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           smcla_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           smcla_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
           
   if (elmset == 'Element PART-1-1.33004'):   # SMCLM TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           smclm_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           smclm_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
          
   if (elmset == 'Element PART-1-1.33006'):   # SMCLP TIB
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           smclp_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           smclp_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
            
  
                      
for i in range(len(time_list)):

    fw.write(" %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n"%(time_list[i],als_CU1_list[i],als_CTF1_list[i],pfl_CU1_list[i],pfl_CTF1_list[i],pcapl_CU1_list[i],pcapl_CTF1_list[i],pcapm_CU1_list[i],pcapm_CTF1_list[i],lcla_CU1_list[i],lcla_CTF1_list[i],lclm_CU1_list[i],lclm_CTF1_list[i],lclp_CU1_list[i],lclp_CTF1_list[i],smcla_CU1_list[i],smcla_CTF1_list[i],smclm_CU1_list[i],smclm_CTF1_list[i],smclp_CU1_list[i],smclp_CTF1_list[i]))
#    fw.write(" %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f   \n"%(time_list[i],als_CU1_list[i],als_CTF1_list[i],pfl_CU1_list[i],pfl_CTF1_list[i],pcapl_CU1_list[i],pcapl_CTF1_list[i],pcapm_CU1_list[i],pcapm_CTF1_list[i]))

odb.close()
