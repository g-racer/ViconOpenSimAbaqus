from odbAccess import *

class OdbReader:
  def __init__(self, odbFile):
    self.odbFile = odbFile
    self.odb = openOdb(path=self.odbFile)
    self.count = 0
    self.copPat1 = []
    self.copPat2 = []
    self.copPat3 = []
    self.copMed1 = []
    self.copMed2 = []
    self.copMed3 = []
    self.copLat1 = []
    self.copLat2 = []
    self.copLat3 = []
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
    fw=open("COP.results",'w')
    fw.write("     TIME    MED_XN1 MED_XN2 MED_XN3  LAT_XN1 LAT_XN2 LAT_XN3  PAT_XN1 PAT_XN2 PAT_XN3      \n")

    numData = len(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN1      SURF-BONE-T_TRAY_M/SURF-BONE-FEM_COMP_M'].data)
    for i in range(numData):
            self.timeList.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN1      SURF-BONE-T_TRAY_M/SURF-BONE-FEM_COMP_M'].data[i][0])
            self.copMed1.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN1      SURF-BONE-T_TRAY_M/SURF-BONE-FEM_COMP_M'].data[i][1])
            self.copMed2.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN2      SURF-BONE-T_TRAY_M/SURF-BONE-FEM_COMP_M'].data[i][1])
            self.copMed3.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN3      SURF-BONE-T_TRAY_M/SURF-BONE-FEM_COMP_M'].data[i][1])
            self.copLat1.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN1      SURF-BONE-T_TRAY_L/SURF-BONE-FEM_COMP_L'].data[i][1])
            self.copLat2.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN2      SURF-BONE-T_TRAY_L/SURF-BONE-FEM_COMP_L'].data[i][1])
            self.copLat3.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN3      SURF-BONE-T_TRAY_L/SURF-BONE-FEM_COMP_L'].data[i][1])
            self.copPat1.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN1      SURF-BONE-DOME/SURF-BONE-FEM_COMP_P'].data[i][1])
            self.copPat2.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN2      SURF-BONE-DOME/SURF-BONE-FEM_COMP_P'].data[i][1])
            self.copPat3.append(stepName.historyRegions['ElementSet  PIBATCH'].historyOutputs['XN3      SURF-BONE-DOME/SURF-BONE-FEM_COMP_P'].data[i][1])
            fw.write(" %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n"%(self.timeList[i],self.copMed1[i],self.copMed2[i],self.copMed3[i],self.copLat1[i],self.copLat2[i],self.copLat3[i],self.copPat1[i],self.copPat2[i],self.copPat3[i]))

    self.odb.close()
if __name__ == '__main__':

  import sys

  file = sys.argv[1]

  OdbReader(file).run()

