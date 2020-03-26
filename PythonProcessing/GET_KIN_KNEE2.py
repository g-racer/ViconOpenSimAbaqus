from odbAccess import *
import sys

# script for extracting KKS actuator data

s=sys.argv[1]

odb = openOdb(path='%s'%(s))

time_list = []
tf_fe_list = []
tf_vv_list = []
tf_ie_list = []
tf_ml_list = []
tf_ap_list = []
tf_si_list = []


fw=open("COMP_KIN_KNEE.results",'w')
fw.write("TIME(sec)    tf_FE    tf_VV    tf_IE   tf_ML   tf_AP   tf_SI   \n")

for elmset in odb.steps['Step-3'].historyRegions.keys():

############################### QUADRICEPS RESULTS ################################
   # quad loads
   if (elmset == 'Element PART-1-1.152'):   # PATELLA LOADS AXIS
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data)
       for i in range(numData):
           time_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][0])
           tf_fe_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CUR1'].data[i][1])
           tf_ml_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])


   if (elmset == 'Element PART-1-1.151'):   # PATELLA LOADS AXIS
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data)
       for i in range(numData):
           tf_vv_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CUR1'].data[i][1])
           tf_ap_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])   

   if (elmset == 'Element PART-1-1.150'):   # PATELLA LOADS AXIS
       numData = len(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CTF1'].data)
       for i in range(numData):
           tf_ie_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CUR1'].data[i][1])
           tf_si_list.append(odb.steps['Step-3'].historyRegions[elmset].historyOutputs['CU1'].data[i][1])   
                      

for i in range(len(time_list)):

    fw.write(" %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  \n"%(time_list[i],tf_fe_list[i],tf_vv_list[i],tf_ie_list[i],tf_ml_list[i],tf_ap_list[i],tf_si_list[i]))

odb.close()
