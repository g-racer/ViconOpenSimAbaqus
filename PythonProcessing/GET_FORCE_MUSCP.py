from odbAccess import *
import sys

# script for extracting KKS actuator data

s=sys.argv[1]

odb = openOdb(path='%s'%(s))

time_list = []


RF_CU1_list = []
RF_CTF1_list = []
VI_CU1_list = []
VI_CTF1_list = []
VL_CU1_list = []
VL_CTF1_list = []
VM_CU1_list = []
VM_CTF1_list = []

GASLAT_CU1_list = []
GASLAT_CTF1_list = []







fw=open("FORCE_MUSCP.results",'w')
fw.write("TIME(sec) recfem_excursion  recfem_force vasint_excursion vasint_force  vaslat_excursion  vaslat_force  vasmed_excursion  vasmed_force  \n")

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


        
                      
for i in range(len(time_list)):

    fw.write(" %9.3f  %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  \n"%(time_list[i],RF_CU1_list[i],RF_CTF1_list[i],VI_CU1_list[i],VI_CTF1_list[i],VL_CU1_list[i],VL_CTF1_list[i],VM_CU1_list[i],VM_CTF1_list[i]))

odb.close()
