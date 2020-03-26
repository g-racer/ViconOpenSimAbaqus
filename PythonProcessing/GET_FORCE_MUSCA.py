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





fw=open("FORCE_MUSCA.results",'w')
fw.write("TIME(sec) gaslat_excursion  gaslat_force gasmed_excursion  gasmed_force  bflh_excursion  bflh_force  bfsh_excursion  bfsh_force   grac_excursion grac_force  semimem_excursion semimem_force semiten_excursion semiten_force  tfl_excursion  tfl_force  sart_excursion sart_force \n")

for elmset in odb.steps['Step-3'].historyRegions.keys():

############################### QUADRICEPS RESULTS ################################
   # quad loads



   if (elmset == 'Element PART-1-1.103900'):   # musc-smemb1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           time_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][0])
           SMEMB_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           SMEMB_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
           
   if (elmset == 'Element PART-1-1.104000'):   # musc-sten1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           STEND_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           STEND_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
          
   if (elmset == 'Element PART-1-1.101000'):   # musc-bflh1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           BFLH_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           BFLH_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])

   if (elmset == 'Element PART-1-1.101100'):   # musc-bfsh1
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           BFSH_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           BFSH_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   	   
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
		   
   if (elmset == 'Element PART-1-1.104200'):   # musc-tfl
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           TFL_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           TFL_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   
   if (elmset == 'Element PART-1-1.101700'):   # musc-gasmed
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           GASMED_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           GASMED_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   
   if (elmset == 'Element PART-1-1.101600'):   # musc-gaslat
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data)
       print numData
       for i in range(numData):
           GASLAT_CU1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])
           GASLAT_CTF1_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data[i][1])
		   

        
        
                      
for i in range(len(time_list)):

    fw.write(" %9.3f  %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f   \n"%(time_list[i],GASLAT_CU1_list[i],GASLAT_CTF1_list[i],GASMED_CU1_list[i],GASMED_CTF1_list[i],BFLH_CU1_list[i],BFLH_CTF1_list[i],BFSH_CU1_list[i],BFSH_CTF1_list[i],GRAC_CU1_list[i],GRAC_CTF1_list[i],SMEMB_CU1_list[i],SMEMB_CTF1_list[i],STEND_CU1_list[i],STEND_CTF1_list[i],TFL_CU1_list[i],TFL_CTF1_list[i],SART_CU1_list[i],SART_CTF1_list[i]))

odb.close()
