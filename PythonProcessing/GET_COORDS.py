# ptest.py

from odbAccess import *

class OdbReader:
  def __init__(self, odbFile):
    self.odbFile = odbFile
    self.odb = openOdb(path=self.odbFile)
    self.count = 0
    self.timeListCFN = []
    self.cfn1List = []
    self.cfn2List = []
    self.cfn3List = []
    self.careaList1 = []
    self.careaList2 = []
    self.cpressList1 = []
    self.cpressList2 = []
    self.LCL = []
    self.SMCL = []
    self.SMCL_OB = []
    self.DMCL = []
    self.PCAPM = []
    self.PCAPL = []
    self.POP = []
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

    numframes=len(self.odb.steps['Step-3'].frames)
    fw=open("COORDS.results",'w')
    fw.write("     TIME    FEM_1x  FEM_1y  FEM_1z  FEM_2x  FEM_2y  FEM_2z  FEM_3x  FEM_3y  FEM_3z  FEM_4x  FEM_4y  FEM_4z  TIB_1x  TIB_1y  TIB_1z  TIB_2x  TIB_2y  TIB_2z  TIB_3x  TIB_3y  TIB_3z  TIB_4x  TIB_4y  TIB_4z  PAT_1x  PAT_1y  PAT_1z  PAT_2x  PAT_2y  PAT_2z  PAT_3x  PAT_3y  PAT_3z  PAT_4x  PAT_4y  PAT_4z  \n")
    
    #cycle through all frames    
    for i in range(numframes):
        Frame= self.odb.steps['Step-3'].frames[i]
        refn = self.odb.rootAssembly.instances['PART-1-1'].nodeSets['GSNODES']
        coordinate = Frame.fieldOutputs['COORD']
        refcoordinate = coordinate.getSubset(region=refn)

        # array of each node X,Y,Z value for the first 4 nodes in GS_SET (tibial nodes)
        for m in refcoordinate.values:
            #print m.nodeLabel,Frame.frameValue,m.data[0],m.data[1],m.data[2]
            self.coord.append(m.data[0]) 
            self.coord.append(m.data[1]) 
            self.coord.append(m.data[2]) 
        
        fw.write(" %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n"%(Frame.frameValue, self.coord[0], self.coord[1], self.coord[2], self.coord[3], self.coord[4], self.coord[5], self.coord[6], self.coord[7], self.coord[8], self.coord[9],self.coord[10], self.coord[11], self.coord[12], self.coord[13], self.coord[14], self.coord[15], self.coord[16], self.coord[17], self.coord[18], self.coord[19],self.coord[20], self.coord[21], self.coord[22], self.coord[23], self.coord[24], self.coord[25], self.coord[26], self.coord[27], self.coord[28], self.coord[29],self.coord[30], self.coord[31], self.coord[32]))
        
        self.coord=[]
        
        
    self.odb.close()


if __name__ == '__main__':

  import sys

  file = sys.argv[1]

  OdbReader(file).run()
