import csv
import os

def parseOutFile(outFile, moleculeName, csdCode):
  # print("in parseOutput File")
  iF = open(outFile, "r")
  oF = open("./summary.csv", "a")
  writer = csv.writer(oF)
  data = [moleculeName,csdCode,'None','None','None','None','None','None','None','None','None','None','None','None']
  fileNameIndex = 0
  csdIndex = 1
  homoSvpIndex = 2
  homoTzvpdIndex = 3
  homoSinglePointSvpIndex = 4
  homoSinglePointTzvpdIndex = 5
  hetroPlusSvpIndex = 6
  hetroPlusTzvpdIndex = 7
  hetroPlusSinglePointSvpIndex = 8
  hetroPlusSinglePointTzvpdIndex = 9
  hetroNegSvpIndex = 10
  hetroNegTzvpdIndex = 11
  hetroNegSinglePointSvpIndex = 12
  hetroNegSinglePointTzvpdIndex = 13
  for line in iF:
    if "CSD code:" in line:
      data[csdIndex] = str(line.split(":")[1].split()[0])
    if "Molecule name:" in line:
      data[fileNameIndex] = str(line.split(":")[1].split()[0])
    if "Bond Energies: (M-Homolysis)" in line:
      # print("found Bond Energies: (M-Homolysis)")
      line = iF.readline()
      line = iF.readline()
      def2SVPLine = line
      def2SVPEnergy = def2SVPLine.split('=')[1].split()[0]
      line = iF.readline()
      def2TZVPDLine = line
      def2TZVPDEnergy = def2TZVPDLine.split('=')[1].split()[0]
      #single point energies
      line = iF.readline()
      spDef2SVPLine = line
      line = iF.readline()
      spDef2TZVPDLine = line
      spDef2SVPEnergy = spDef2SVPLine.split('=')[1].split()[0]
      spDef2TZVPDEnergy =  spDef2TZVPDLine.split('=')[1].split()[0]
      #write to the array
      data[homoSvpIndex] = str(def2SVPEnergy)
      data[homoTzvpdIndex] = str(def2TZVPDEnergy)
      data[homoSinglePointSvpIndex] = str(spDef2SVPEnergy)
      data[homoSinglePointTzvpdIndex] = str(spDef2TZVPDEnergy)
    if "Bond Energies: ((M+) Heterolysis)" in line:
      line = iF.readline()
      line = iF.readline()
      def2SVPLine = line
      def2SVPEnergy = def2SVPLine.split('=')[1].split()[0]
      line = iF.readline()
      def2TZVPDLine = line
      def2TZVPDEnergy = def2TZVPDLine.split('=')[1].split()[0]
      #single point energies
      line = iF.readline()
      spDef2SVPLine = line
      line = iF.readline()
      spDef2TZVPDLine = line
      spDef2SVPEnergy = spDef2SVPLine.split('=')[1].split()[0]
      spDef2TZVPDEnergy =  spDef2TZVPDLine.split('=')[1].split()[0]
      #write to array
      data[hetroPlusSvpIndex] = str(def2SVPEnergy)
      data[hetroPlusTzvpdIndex] = str(def2TZVPDEnergy)
      data[hetroPlusSinglePointSvpIndex] = str(spDef2SVPEnergy)
      data[hetroPlusSinglePointTzvpdIndex] = str(spDef2TZVPDEnergy)
    if "Bond Energies: ((M-) Heterolysis)" in line:
      line = iF.readline()
      line = iF.readline()
      def2SVPLine = line
      def2SVPEnergy = def2SVPLine.split('=')[1].split()[0]
      line = iF.readline()
      def2TZVPDLine = line
      def2TZVPDEnergy = def2TZVPDLine.split('=')[1].split()[0]
      #single point energies
      line = iF.readline()
      spDef2SVPLine = line
      line = iF.readline()
      spDef2TZVPDLine = line
      spDef2SVPEnergy = spDef2SVPLine.split('=')[1].split()[0]
      spDef2TZVPDEnergy =  spDef2TZVPDLine.split('=')[1].split()[0]
      #write to array
      data[hetroNegSvpIndex] = str(def2SVPEnergy)
      data[hetroNegTzvpdIndex] = str(def2TZVPDEnergy)
      data[hetroNegSinglePointSvpIndex] = str(spDef2SVPEnergy)
      data[hetroNegSinglePointTzvpdIndex] = str(spDef2TZVPDEnergy)
  # print("data before writing: " + str(data))
  writer.writerow(data)

def startCsvFile():
  # print("in csv file start")
  oF = open("./summary.csv", "w")
  writer = csv.writer(oF)
  header = ['File Name','CSD code','Homo Def2-SVP','Homo Def2-TZVPD', 'Homo-Single-Point Def2-SVP', 'Homo-Single-Point Def2-TZVPD',
    'Hetero+ Def2-SVP','Hetero+ Def2-TZVPD', 'Hetero+ Single-Point Def2-SVP', 'Hetero+ Single-Point Def2-TZVPD',
    'Hetero- Def2-SVP', 'Hetero- Def2-TZVPD', 'Hetero- Single-Point Def2-SVP', 'Hetero- Single-Point Def2-TZVPD']
  # print("header: " + str(header))
  writer.writerow(header)
  oF.close()

def getMoleculeInfo(folder):
  # print("folder: " + str(folder))
  count = 0 
  for file in os.listdir(folder):
    print(file + "\n")
    if "." in file:
      if file.split(".")[1] == "xyz":
        print("found the xyz file")
        moleculeName = folder.split(".")[0].split()[0]
        iF = open(folder + "/" + file)
        iF.readline()
        line = iF.readline()
        csdCode = line.split(";")[3]
        print("molecule " + str(count) + moleculeName + " ,csd: " + csdCode + "\n")
        iF.close()
    count += 1
  return moleculeName, csdCode

# def writeCSV():
#   startCsvFile()
#   path = "./"
#   for folder in os.listdir("./"):
#     if folder == "summary.csv":
#       continue
#     path = "./" + folder
#     for sub_folder in os.listdir("./" + folder + "/"):
#       path = "./" + folder + "/" + sub_folder + "/"
#       if sub_folder == "completed":
#         path = "./" + folder + "/" + sub_folder + "/"
#         for file in os.listdir(path):
#           if file.split(".")[1] == "out":
#             parseOutFile(path + file)