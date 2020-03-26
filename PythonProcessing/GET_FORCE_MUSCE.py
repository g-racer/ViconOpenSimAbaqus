from odbAccess import *
import sys

# script for extracting KKS actuator data

s=sys.argv[1]

odb = openOdb(path='%s'%(s))

time_list = []
GASLAT_CU1_list = []
GASLAT_CTF1_list = []
GASMED_CU1_list = []
GASMED_CTF1_list = []


RF_CU1_list = []
RF_CTF1_list = []
VI_CU1_list = []
VI_CTF1_list = []
VL_CU1_list = []
VL_CTF1_list = []
VM_CU1_list = []
VM_CTF1_list = []

BFLH_CU1_list = []
BFLH_CTF1_list = []
BFSH_CU1_list = []
BFSH_CTF1_list = []

GRAC_CU1_list = []
GRAC_CTF1_list = []
SMEMB_CU1_list = []
SMEMB_CTF1_list = []
STEND_CU1_list = []
STEND_CTF1_list = []

TFL_CU1_list = []
TFL_CTF1_list = []

SART_CU1_list = []
SART_CTF1_list = []





fw=open("FORCE_MUSCE.results",'w')
fw.write("TIME(sec)    gasmed_excursion  gasmed_force recfem_excursion  recfem_force vasint_excursion vasint_force  vaslat_excursion  vaslat_force  vasmed_excursion  vasmed_force    grac_excursion grac_force  semimem_excursion semimem_force semiten_excursion semiten_force  sart_excursion sart_force \n")

for elmset in odb.steps['Step-3'].historyRegions.keys():

############################### QUADRICEPS RESULTS ################################
   # quad loads
   if (elmset == 'Element PART-1-1.1500'):   # musc-recfem1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           time_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][0])
           RF_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           RF_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])

   if (elmset == 'Element PART-1-1.2500'):   # musc-vasmed1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           VM_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           VM_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])

   if (elmset == 'Element PART-1-1.3500'):   # musc-vaslat1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           VL_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           VL_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
                      
   if (elmset == 'Element PART-1-1.4500'):   # musc-vasint1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           VI_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           VI_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])

   if (elmset == 'Element PART-1-1.103900'):   # musc-smemb1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           SMEMB_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           SMEMB_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
           
   if (elmset == 'Element PART-1-1.104000'):   # musc-sten1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           STEND_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           STEND_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
          	   
   if (elmset == 'Element PART-1-1.102800'):   # musc-grac
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           GRAC_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           GRAC_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   
   if (elmset == 'Element PART-1-1.103800'):   # musc-sart
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           SART_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           SART_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   	   
   if (elmset == 'Element PART-1-1.101700'):   # musc-gasmed
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           GASMED_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           GASMED_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   


        
        
                      
for i in range(len(time_list)):

    fw.write(" %9.3f %9.3f %9.3f  %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n"%(time_list[i],GASMED_CU1_list[i],GASMED_CTF1_list[i],RF_CU1_list[i],RF_CTF1_list[i],VI_CU1_list[i],VI_CTF1_list[i],VL_CU1_list[i],VL_CTF1_list[i],VM_CU1_list[i],VM_CTF1_list[i],GRAC_CU1_list[i],GRAC_CTF1_list[i],SMEMB_CU1_list[i],SMEMB_CTF1_list[i],STEND_CU1_list[i],STEND_CTF1_list[i],SART_CU1_list[i],SART_CTF1_list[i]))

odb.close()
