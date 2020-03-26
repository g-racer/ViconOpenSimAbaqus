# ptest.py

from odbAccess import *

class OdbReader:
  def __init__(self, odbFile):
    self.odbFile = odbFile
    self.odb = openOdb(path=self.odbFile)
    self.count = 0
    self.cpressList1 = []
    self.cpressList2 = []
    self.cpressList3 = []
    self.cpressList4 = []
    self.timeList = []
    self.coord = []

  def run(self):
    def unit(vec):
      unitvec = []
      pre_mag = vec[0]**2 + vec[1]**2 + vec[2]**2
      mag = pre_mag**0.5
      if mag == 0.0:
         mag = 0.000000001
      for i in range (len(vec)):
        unitvec.append(vec[i]/mag)
      return unitvec

    def cross_product(v1, v2):
      x3 = v1[1]*v2[2] - v1[2]*v2[1]
      y3 = v1[2]*v2[0] - v1[0]*v2[2]
      z3 = v1[0]*v2[1] - v1[1]*v2[0]
      v = [x3, y3, z3]
      return v

    def dot_product(d1, d2):
      x = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]
      return x
    
    stepName = self.odb.steps['Step-3']
    count = 1

    fw=open("CPRESS.results",'w')
    fw.write("     TIME    TF_MED_PRESS   TF_LAT_PRESS   PF_PRESS \n")
    
    for frame in stepName.frames:
        count=count + 1
        maxCPress2 = 0
        cpressSet2 = frame.fieldOutputs['CPRESS   SURF-BONE-T_TRAY_M/SURF-BONE-FEM_COMP_M']
        maxCPress3 = 0
        cpressSet3 = frame.fieldOutputs['CPRESS   SURF-BONE-T_TRAY_L/SURF-BONE-FEM_COMP_L']
        maxCPress4 = 0
        cpressSet4 = frame.fieldOutputs['CPRESS   SURF-BONE-DOME/SURF-BONE-FEM_COMP_P']
              
        
        for cpressValue2 in cpressSet2.values:
          if (cpressValue2.data > maxCPress2):
            maxCPress2 = cpressValue2.data
        self.cpressList2.append(maxCPress2)
        
        for cpressValue3 in cpressSet3.values:
          if (cpressValue3.data > maxCPress3):
            maxCPress3 = cpressValue3.data
        self.cpressList3.append(maxCPress3)

        for cpressValue4 in cpressSet4.values:
          if (cpressValue4.data > maxCPress4):
            maxCPress4 = cpressValue4.data
        self.cpressList4.append(maxCPress4)
        

        fw.write(" %9.3f %9.3f %9.3f %9.3f \n"%(frame.frameValue, maxCPress2, maxCPress3, maxCPress4 ))

    self.odb.close()
if __name__ == '__main__':

  import sys

  file = sys.argv[1]

  OdbReader(file).run()
