#!/bin/env python3
#Written By Braden Borough 

import sys, os,re, shutil, numpy as np, time, subprocess
import os.path


#This list is currently incomplete But it is a mapping of metals and nonmetals to their atomic number
metals_lib = {"26":"FE","51":"SB","75":"RE","3":"LI","4":"BE","11":"NA","12":"MG",
	"78":"PT","23":"V","74":"W","27":"CO","29":"CU","46":"PD","44":"RU","30":"ZN",
	"47":"AG","48":"CD","21":"SC","45":"RH","28":"NI","42":"MO","25":"MN","77":"IR",
	"79":"AU","80":"HG","76":"OS","40":"ZR","39":"Y","73":"TA","24":"CR","22":"TI",
	"43":"TC","105":"HA","57":"LA"}
non_metals_lib = {"7":"N","15":"P","1":"H","6":"C","16":"S","8":"O","9":"F","10":"NE",
	"17":"CL","18":"AR","34":"SE","35":"BR","36":"KR","53":"I","54":"XE","86":"RN"}
# ch3_heterolysisPositiveMetal_lib = {"electronic":"-39.7189771277","zero-point":"-39.691062","enthalpy":"-39.687244","free":"-2338.424417"}
# ch3_heterolysisNegativeMetal_lib = {"electronic":"-39.3881869409","zero-point":"-39.356948","enthalpy":"-39.353147","free":"-39.374372"}
# ch3_homolysis_lib = {"electronic":"-39.7532353865","zero-point":"-39.723615","enthalpy":"-39.719575","free":"-39.74176"}
ch3_Def2_SVP_lib = {"homolysis":"-39.7355115751","heterolysisPositiveMetal":"-39.3768993516","heterolysisNegativeMetal":"-39.7051801432"}
ch3_Def2_TZVPD_lib = {"homolysis":"-39.7916491654","heterolysisPositiveMetal":"-39.4305274697","heterolysisNegativeMetal":"-39.7870838290"}
metal_heterolysisPositive_lib = {}
metal_heterolysisNegative_lib = {}
metal_homolysis_lib = {}
metal_MetalAndCh3_lib = {}

#Loads defaults from the default file into the inputs variable
def default():
	inputs = {}
	iF = open(os.path.expanduser("~/TSS/bin/.default"), "r")
	line = iF.readline()
	while line:
		if "spin" in line:
			inputs["spin"] = line.split(':')[1]
		elif "charge" in line:
			inputs["charge"] = line.split(':')[1]
		elif "basis" in line:
			inputs["basis"] = line.split(':')[1]
		elif "method" in line:
			inputs["method"] = line.split(':')[1]
		elif "batch" in line:
			inputs["batch"] = line.split(':')[1]
		elif "metal" in line:
			inputs["mbasis"] = line.split(':')[1]
		elif "denfit" in line:
			inputs["denfit"] = line.split(':')[1]
		#check with Taylor, what does the crest_selectivity do does it need to be removed from default???
		line = iF.readline()
	iF.close()
	inputs["opt"] = "modred"
	return inputs


#Parses the input file and overrides the defaults with values found in the input file
def parseInput(inputFile, inputs, xyz_file):
	moleculeName = xyz_file.split(".")[0]
	name = inputFile.split('.')[0]
	parseInputTestFile = open("parseInputTestFile.txt", "w")
	parseInputTestFile.write("argv[2].split('/')[0] + xyz: " + sys.argv[2].split("/")[-1] + "\n")
	# parseInputTestFile.write("xyz file: " + xyz_file + "\n")
	# parseInputTestFile.write("name from .in file: " + name + "\n")
	# parseInputTestFile.write("name from .xyz file: " + xyz_file.split('.')[1] + "\n")
	extension = inputFile.split('.')[1]
	if extension == "in":
		pass
	else:
		sys.exit("input file wasn't a .in file")
	iF = open(inputFile, "r")
	coords = []
	line = iF.readline()
	while line:
		if "basis" in line:
			inputs["basis"] = line.split(':')[1]
		 #Required - Default in file
		elif "method" in line:
			inputs["method"] = line.split(':')[1]
		 #Required - Gaussian has a default
		elif "temperature" in line:
			inputs["temperature"] = line.split(':')[1]
		#Not Required, Gaussian Defaults to Gas
		elif "solvent" in line:
			inputs["solvent"] = line.split(':')[1]
		elif "metal" in line:
			inputs["mbasis"] = line.split(':')[1]
		elif "solvent_model" in line:
			inputs["solvent_model"] = line.split(':')[1]
		 #Required - Default in file
		elif "denfit" in line:
			inputs["denfit"] = line.split(':')[1]
		elif "memory" in line:
			inputs["mem"] = line.split(':')[1]
		elif "memper" in line:
			inputs["memper"] = line.split(':')[1]
		elif "confs" in line:
			inputs["confs"] = line.split(':')[1]
			os.environ['NUM_CONFS_REQUESTED'] = inputs["confs"]
		elif "extra" in line:
			inputs["extra"] = line.split(':')[1]
		elif "difficulty" in line:
			inputs["difficulty"] = line.split(':')[1]
		elif "conformational_leniency" in line:
			inputs["leniency"] = line.split(':')[1]
		elif "run_crest" in line: #Check with Taylor???
			inputs["run_crest"] = line.split(':')[1]
			os.environ['RUN_CREST'] = inputs["run_crest"] #??? using env variable to
		elif "vibration_test" in line:
			inputs["vibration_test"] = line.split(':')[1]
			os.environ['VIBRATION_TEST'] = inputs["vibration_test"]
			line = iF.readline()
			inputs["atom1"] = str(line.split(':')[1])
			os.environ['ATOM1'] = inputs["atom1"]
			line = iF.readline()
			inputs["atom2"] = str(line.split(':')[1])
			os.environ['ATOM2'] = inputs["atom2"]
		elif "basic_optimization" in line:
			inputs["opt"] = "basicOpt"
			os.environ['BASIC_OPTIMIZATION'] = line.split(':')[1] #shoulc be set to yes if the user only wants a basic optimization job
		elif "run_type" in line:
			os.environ["RUN_TYPE"] = line.split(":")[1] 
		line = iF.readline()
	iF.close()
	#coords = getCoords(inputs,xyz_file)
	return inputs# coords

#Gets the coordinates from the xyz file and puts them into a list
# def getCoords(inputs,xyz_file):
# 	coords = []
# 	#iF = open(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "r")
# 	iF = open(xyz_file,"r")
# 	moleculeName = xyz_file.split(".")[0]
# 	os.environ['MOLECULE_NAME'] = str(moleculeName) #set the name of the molecule in the environment for later use
# 	line = iF.readline()
# 	line = iF.readline()
# 	parseFreezes(line)
# 	line = iF.readline()
# 	while line:
# 		coords.append(line)
# 		line = iF.readline()
# 	return coords
	
def addComBasisSet(file_name, inputs, check_file_name):
	oF = open(file_name, "a")
	oF.write("--link1--\n")
	# oF.write("%nprocshared=7\n")
	# oF.write("%mem=" + str(int(inputs["memper"]) // int(inputs["confs"])) + "GB\n")
	iF = open("../basisSetTemplate.com", "r")
	line = iF.readline()
	while "link1" not in line:
		line = iF.readline()
	info = iF.readlines()
	for line in info:
		if "0 1\n" in line:
			oF.write(inputs["charge"].strip() + " " + inputs["spin"] + "\n")
			continue
		if "%chk=" in line:
			oF.write(check_file_name)
			continue
		oF.write(line)
	oF.close()

#Builds a com file using values from input file and defaults for any not specified
def buildCom(inputs, coords, f_name):
	# oF = open(f_name, 'w')
	# oF.write("%nprocshared=7\n")
	# oF.write("%mem="+ str(int(inputs["memper"])//int(inputs["confs"])) + "GB\n") #insert into middle of this if fails: str(int(inputs["memper"])//int(inputs["confs"]))
	# basic_opt = os.getenv('BASIC_OPTIMIZATION')
	single_point_calculations = os.getenv('SINGLE_POINT')
	if single_point_calculations != None: #change back to basic_opt if failure
		if "yes" in single_point_calculations:
			com_file_name = f_name.split(".")[0] + "-singlePoint.com"
			oF = open(f_name.split(".")[0] + "-singlePoint.com", 'w')
			oF.write("%nprocshared=28\n")
			total_mem = int(inputs["memper"]) - 2
			oF.write("%mem="+ str(total_mem) + "GB\n")
			# oF.write("%mem="+ str(int(inputs["memper"])//int(inputs["confs"])) + "GB\n") #insert into middle of this if fails: str(int(inputs["memper"])//int(inputs["confs"]))
			basic_opt = os.getenv('BASIC_OPTIMIZATION')
			# oF.write("#opt " + inputs["method"].strip() + "/" + "genecp"  + " freq=noraman" + " integral=ultrafine" + " " + inputs["extra"].strip()) #insert after opt if fails: " freq=noraman " + 
			# oF.write("#" + inputs["method"].strip() + "/" + "genecp" + " integral=ultrafine") #+ " " + inputs["extra"].strip()) #insert after opt if fails: " freq=noraman " +
			check_file_name = "%chk=" + f_name.split(".")[0] + ".chk\n"
			oF.write(check_file_name) # oF.write("%chk=checkFile.chk\n") 
			oF.write("#mn15" + "/" + "genecp" + " integral=ultrafine") #+ " " + inputs["extra"].strip()) #insert after opt if fails: " freq=noraman " + 
	else:
		# oF.write("#opt=" + inputs["opt"].strip() + " freq=noraman " + inputs["method"].strip() + "/" + "genecp"  + " integral=ultrafine" + " " + inputs["extra"].strip())
		com_file_name = f_name
		oF = open(f_name, 'w')
		oF.write("%nprocshared=16\n")
		total_mem = int(inputs["memper"]) - 2
		oF.write("%mem="+ str(total_mem) + "GB\n")
		# oF.write("%mem="+ str(int(inputs["memper"])//int(inputs["confs"])) + "GB\n") #insert into middle of this if fails: str(int(inputs["memper"])//int(inputs["confs"]))
		check_file_name = "%chk=" + f_name.split(".")[0] + ".chk\n"
		oF.write(check_file_name) # oF.write("%chk=checkFile.chk\n")
		basic_opt = os.getenv('BASIC_OPTIMIZATION')
		oF.write("#opt" + " mn15" + "/" + "genecp"  + " integral=ultrafine") #+ " " + inputs["extra"].strip()) #insert after opt if fails: " freq=noraman " + 
	oF.write("\n\n")
	# oF.write("TSS")
	csd_Code = os.environ.get('CSD_CODE')
	oF.write(str(csd_Code))
	oF.write("\n\n")
	oF.write(inputs["charge"].strip() + " " + inputs["spin"] + "\n")
	for coord in coords:
		coord = coord.strip()
		if coord[0].isalpha():
			oF.write(coord + "\n")
		else:
			coord = coord.split()
			if coord[0] in metals_lib:
				string = " "
				coord[0] = metals_lib[str(coord[0])]
				oF.write(string.join(coord)+ "\n")
			else:
				string = " "
				coord[0] = non_metals_lib[str(coord[0])]
				oF.write(string.join(coord)+ "\n")
	oF.write("\n")
	writeFreezes(oF, coords, inputs)
	oF.write("\n")
	oF.write("\n")
	oF.close()
	addComBasisSet(com_file_name, inputs, check_file_name)

#Writes the freezes in the modred file
def writeFreezes(outFile,coords, inputs):
	#Freezes = []
	#libFile = open(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "r")
	#libFile.readline()
	#freeze_line = libFile.readline()
	#freeze_line = freeze_line[2:]
	#freeze_line = freeze_line.split(';')
	#freeze_line.pop(-1)
	#if "modred" in inputs["opt"]:
        #	for i in range(0, len(freeze_line)):
        #        	val = freeze_line[i].split('-')
        #        	Freezes.append(val[0])
        #        	Freezes.append(val[1])
        #	outFile.write("\n")
	#for i in range(0,len(Freezes)-1,2):
        #	outFile.write("B " + str(int(Freezes[i]) + 1) + " " + str(int(Freezes[i+1]) + 1) + " F\n")
	basic_opt = os.getenv("BASIC_OPTIMIZATION")
	if "modred" in inputs["opt"]:
		if basic_opt != None:
			if "yes" not in basic_opt:
				outFile.write("\n")
				for bond in inputs["bonds"]:
					outFile.write("B " + str(bond[0]) + " " + str(bond[1]) + " F\n")
			# else:
			# 	outFile.write("\n")
		else:
			outFile.write("\n")
			for bond in inputs["bonds"]:
				outFile.write("B " + str(bond[0]) + " " + str(bond[1]) + " F\n")
	writeGenecp(outFile, coords, inputs)	


#Writes the genecp section in the modred and TS search file
def writeGenecp(outFile,coords, inputs):
	metals, non_metals = getAtomTypes(coords)
	if "modred" in inputs["opt"]:
		outFile.write("\n")
	for val in metals:
		outFile.write(val + " ")
	outFile.write("0\n")
	outFile.write(inputs["mbasis"].strip() + "\n****\n")
	for val in non_metals:
		outFile.write(val + " ")
	outFile.write("0\n" + inputs["basis"].strip() + "\n****\n\n")
	for val in metals:
		outFile.write(val + " ")
	outFile.write("0\n")
	outFile.write(inputs["mbasis"].strip())
		


#Gets a list of what types of atoms are in the file to be used in the write genecp section
def getAtomTypes(coords):
	metals = set()
	non_metals = set() 
	for coord in coords:
		coord = coord.strip()
		if coord.split()[0].upper() in metals_lib.values():
			metals.add(coord.split()[0])
		elif coord.split()[0] in metals_lib:
			metals.add(metals_lib[str(coord.split()[0])])
		else:
			if coord[0].isalpha():
				non_metals.add(coord.split()[0])
			else:
				non_metals.add(non_metals_lib[str(coord.split()[0])])
	return metals, non_metals
	

#Copies the base_input as defined in the input file and modifies it to build a new xyz file
#Currently only does subractions
#Currently not being used
def buildLibraryInputs(lib_location):
	shutil.copy(os.path.expanduser("~/TSS/libs/base_templates/" + inputs["library"].strip()), "temp.xyz")
	tempF = open("temp.xyz", "r")
	finalF = open("crest_conformers.xyz", "w")
	atomNumber = int(tempF.readline())
	atomNumber -= len(inputs["subtract"])
	#need to include adding and substituting
	tempF.readline()	
	finalF.write(str(atomNumber) + "\n\n")
	for i in range(1, atomNumber + 1):
		line = tempF.readline()
		if str(i) in inputs["subtract"]:
			pass
		#elif i in inputs["substitute"]:
		#elif i in inputs["add"]:
		else:
			finalF.write(line)

#Uses the results of the modred to make com files for the modred section
def modredCrest(crest_file, inputs):
	# print("running modredCrest")
	energies = []
	acceptable_energy = True
	os.chdir("modred")
	coords_list = []
	names_list = []
	num_structures = 0
	if os.path.isfile("../"+crest_file):
		iF = open("../"+crest_file,"r")
	elif os.path.isfile(crest_file):
		iF = open(crest_file, "r")
	#iF = open("../"+crest_file,"r")
	line = iF.readline()
	while line:
		line = iF.readline()
		energies.append(float(line)*627.51)
		line = iF.readline()
		coords = []
		while line:
			if line.strip()[0].isalpha():
				coords.append(line)
				line = iF.readline()
			else:
				break
		#for i in range(0,len(energies)-1):
		#	if energies[-1] -energies[i] < 0.10:
		#		acceptable_energy = False
		#		energies.pop()
		#		break
		if acceptable_energy:
			coords_list.append(coords)
			names_list.append("conf" + str(num_structures) + ".com")
			num_structures +=1
		acceptable_energy = True
	inputs["numconfs"] = num_structures
	#for i in range(0,len(coords_list)):
#		buildCom(inputs, coords_list[i], names_list[i])
	try:
		for i in range(len(coords_list)-1, len(coords_list)-int(inputs["confs"])-1,-1):
			#continue
			buildCom(inputs, coords_list[i], names_list[i]) #??? not sure about commenting this out
	except:
		print("okay\n")
	os.chdir("../")
# def modredRangeCreation():
	#This will take the current modreds and create multiple ones with differing frozen bond lengths

#this is the modredCrest function but altered to accept only the original xyz file. This is when the user
#does not wish to run crest
def modredNoCrest(xyz_file, inputs):
	#energies = []
	acceptable_energy = True
	os.chdir("modred")
	coords_list = []
	names_list = []
	num_structures = 0
	if os.path.isfile("../"+xyz_file):
		iF = open("../"+xyz_file,"r")
	elif os.path.isfile(xyz_file):
		iF = open(xyz_file, "r")
	#iF = open("../"+crest_file,"r")
	line = iF.readline()
	while line:
		line = iF.readline()
		#energies.append(float(line)*627.51)
		line = iF.readline()
		coords = []
		while line:
			if line.strip()[0].isalpha():
				coords.append(line)
				line = iF.readline()
			else:
				break
		#for i in range(0,len(energies)-1):
		#	if energies[-1] -energies[i] < 0.10:
		#		acceptable_energy = False
		#		energies.pop()
		#		break
		if acceptable_energy:
			file_name = str((xyz_file.split("/")[-1]).split(".")[0])
			coords_list.append(coords)
			names_list.append("conf" + file_name + ".com") #str(num_structures)
			num_structures +=1
		acceptable_energy = True
	inputs["numconfs"] = num_structures
	#for i in range(0,len(coords_list)):
#		buildCom(inputs, coords_list[i], names_list[i])
	try:
		for i in range(len(coords_list)-1, len(coords_list)-int(inputs["confs"])-1,-1):
			#continue
			buildCom(inputs, coords_list[i], names_list[i]) #??? not sure about commenting this out
	except:
		print("okay\n")
	os.chdir("../")
# def modredRangeCreation():
	#This will take the current modreds and create multiple ones with differing frozen bond lengths

#Takes the log of the modred and gets the resulting xyz coords
def logtoxyz(f_name):
	inFile = open(f_name, 'r')
	iF = inFile.readlines()
	myLine = 0
	for i, line in enumerate(iF):
		if 'Standard orientation' in line:
			myLine = i
	coords = []
	done = False
	i = myLine + 5
	myRegex = r'\s*\d*\s*(\d*)\s*\d*\s*(.*\s*.*\s*.*)'
	while not done:
		if '--' in iF[i]:
			break
		l = re.findall(myRegex, iF[i], flags=0)	
		line = str(l[0][0]) + '\t' + str(l[0][1])
		coords.append(line)
		i += 1
	inFile.close()
	z = 0
	return coords


#Runs the gaussian jobs and monitors them. If the modred finishes but doesn't have an imaginary vibration it is killed otherwise it runs the transition state search
def gaussianProcesses(inputs):
	# print("in gaussianProcesses")
	commands = []
	switched = []
	optType = []	
	processes = []
	file_names = []
	args = sys.argv
	allDone = False
	os.chdir("modred")
	crest_file = open("../mr_test_file.txt", 'w')
	basic_opt = os.getenv("BASIC_OPTIMIZATION")
	for file in os.listdir(os.getcwd()):
		if file.split(".")[1] != "com": #checks to make sure we are only adding new commands to the array if it's a com file
			continue  
		crest_file.write("in gaussianProcesses: line 429...\n")
		# print("file name to add: " + str(file.split('.')[0]))
		file_names.append(str(file.split('.')[0]))
		# print("filename" + file)
		crest_file.write("filename" + file + "/n")
		commands.append(['/apps/gaussian16/B.01/AVX2/g16/g16', file])
		if basic_opt != None:
			if "yes" in basic_opt:
				optType.append("basicOpt")
			else:
				optType.append("modred")
		else:
			optType.append("modred")
		switched.append(0)
	crest_file.write("line 441...\n")
	basicOptCount = 0
	for com in commands:
		if basic_opt != None:
			if "yes" in basic_opt:
				#os.chdir("../basicOpt/")
				#command = "../modred" + com
				processes.append(subprocess.Popen(com, stdout = subprocess.PIPE, stdin = subprocess.PIPE))
				# switched[basicOptCount] = 1
				# basicOptCount += 1
		else:
			processes.append(subprocess.Popen(com, stdout = subprocess.PIPE, stdin = subprocess.PIPE))
	crest_file.write("line 453...\n")
	while not allDone:
		#if we are only doing a basic optimization we want to skip to the end of this function and not run another gauussian job
		basic_opt = os.getenv("BASIC_OPTIMIZATION")
		allDone = True
		time.sleep(300)
		i = 0
		drawStatus(file_names, processes, optType, switched)
		#print("testing new line")
		run_crest = os.getenv('RUN_CREST')
		# if basic_opt != None:
		# 	if "yes" in basic_opt:
				# os.chdir("../")
				#break
		crest_file.write("line 467...\n")
		for p in processes:
			if p.poll() is None:
				allDone = False
			else:
				if (not switched[i]):
					#hasNeg = checkNegVib(file_names[i] + ".log")
					completed = checkCompleted(file_names[i] + ".log")[0] #!!!CHANGED
					#if hasNeg:
					if completed:
						allDone = False
						coords = logtoxyz(file_names[i] + ".log")
						# os.chdir('../gaussianTS/')
						# if run_crest == "yes\n":
						# 	buildCom(inputs, coords, file_names[i] + ".com") #check with Taylor ??? original code was just this line, no if/elif block
						# elif i == 0:
						# 	buildCom(inputs, coords, file_names[i] + ".com") #check with Taylor ???
						inputs["numconfs"] = findAliveProcesses(processes)
						# processes[i] = subprocess.Popen(['/apps/gaussian16/B.01/AVX2/g16/g16', (file_names[i] + ".com")])
						#optType[i] = "TS Calc"
						optType[i] = "basicOpt"
						# os.chdir('../modred')
						switched[i] = 1
					else:
						optType[i] = "killed"
			i += 1
			crest_file.write("line 493...\n")
	drawStatus(file_names, processes, optType, switched)
	crest_file.write("line 495...\n")
	files = [f for f in os.listdir('.') if f.split('.')[1] is "log"]
	crest_file.write("line 506...\n")
	# print("all done")

#Helper function to get the number of processes still alive
def findAliveProcesses(processes):
	num = 0
	for p in processes:
		if p.poll() is None:
			num +=1
	if num == 1:
		return 2
	return num

#Helper function to make all the directories
def makeDirectories():
	os.mkdir("modred")
	# os.mkdir("gaussianTS")
	# os.mkdir("crest")
	os.mkdir("completed")
	os.mkdir("newMetalXYZs")

#Helper function to make and update the status file to see the status of the jobs
def drawStatus(file_names, processes, optType, switched):
	os.chdir('../')
	statusFile = open("status", "w")
	i = 0
	for p in processes:
		if p.poll() is None:
				statusFile.write(str(file_names[i]) + "          ------> Running " + str(optType[i]) + "\n\n")
		else:
				#Check for negative frequencyi
			if(not switched[i]):
					if optType[i] is "killed":
							statusFile.write(str(file_names[i]) + "          ------> Killed - Log file does not have normal termination.\n\n")
					else:
							statusFile.write(str(file_names[i]) + "          ------> Transitioning from modred to TS Calc\n\n")
			else:
					statusFile.write(str(file_names[i]) + "          ------> Done\n\n")
		i += 1 
	statusFile.close()
	os.chdir('modred')

#Helper function to check the files for imaginary vibrations
def checkNegVib(inFile):
	iF = open(inFile, "r")
	line = iF.readline()
	while line:
		if "Frequencies --" in line:
			Freq = float(line.split()[2])
			if Freq < 0:
				return float(line.split() [3]) > 0
		line = iF.readline()
	iF.close()
	return False
#Helper function to check to see if the file finished without errors	
def checkCompleted(inFile):
	#!!!!!!!Todo: check for an error line, if you find one, don't restart, otherwise, if it gets to the bo=ttom of this function, restart the job
	iF = open(inFile, "r")
	singlePoint = False
	restart = False
	if "singlePoint" in inFile:
		singlePoint = True
	if not singlePoint:
		for line in iF:
			if "Stationary point found" in line:
				return [True, restart] #!!!CHANGED
			elif "Error" in line: #restart the job if it terminated in an error
				restart = True
			pass
	else:
		for line in iF:
			if "SCF Done" in line:
				return [True, restart] #!!!CHANGED
			elif "Error" in line: #Restart the job if it ternminated in an error
				restart = True
			pass
	last_line = line
	iF.close()
	retArr = []
	retArr.append(False, restart) #return an array where 0 is the success, and 1 is the restart bool
	# if "Normal termination" in last_line:
	return retArr#False


#Helper function to make the output for each file
def outputFunc(f_name):
	thval = ""
	iF = open(f_name, 'r')
	line = iF.readline()
	while line:
		if "Zero-point correction" in line:
			no_of_lines = 7
			lines = line
			for i in range(no_of_lines):
				lines+=iF.readline()
			thval = lines
		line = iF.readline()
	iF.close()
	return thval

#moves the .log files that were sucessfull from the modred folder to the completed folder
def moveLogFiles():
	#os.chdir("../")
	successFullRunTypes = []
	source_path = "modred/"
	target_path = "completed/"
	for filename in os.listdir("modred"):
		if filename.split('.')[1] == "log":
			completed = checkCompleted(source_path + filename)
			if completed[0]: #!!!CHANGED
				singlePointBool = False
				shutil.move(source_path + filename, target_path + filename)
				if "singlePoint" in filename:
					singlePointBool = True
				if "homolysis" in filename:
					if singlePointBool:
						successFullRunTypes.append("homolysis-singlePoint")
					else:
						successFullRunTypes.append("homolysis")
				if "heterolysisPositive" in filename:
					if singlePointBool:
						successFullRunTypes.append("heterolysisPositive-singlePoint")
					else:
						successFullRunTypes.append("heterolysisPositive")
				if "heterolysisNegative" in filename:
					if singlePointBool:
						successFullRunTypes.append("heterolysisNegative-singlePoint")
					else:
						successFullRunTypes.append("heterolysisNegative")
			# elif completed[1]:
			# 	#!!!this is where we will restart the job
			# 	pass
	return successFullRunTypes

#fills out the important bond energy summary info on the final output file called "bradenbond.out"
#fills out the important bond energy summary info on the final output file called "bradenbond.out"
def writeOutputSummaryFile(output_file, sucessfulRunTypes):
	run_type = os.getenv("RUN_TYPE")
	keyWords = []
	if "heterolysisPositive" in run_type:
		keyWords.append("heterolysisPositive")
	if "heterolysisNegative" in run_type:
		keyWords.append("heterolysisNegative")
	if "homolysis" in run_type:
		keyWords.append("homolysis")
	
	oF = open(output_file, "a")
	csd_Code = os.environ.get('CSD_CODE')
	oF.write("CSD code: " + str(csd_Code) + "\n" + "Molecule name: " + str(os.environ.get("MOLECULE_NAME")) + "\n")
	if "homolysis" in keyWords:
		oF.write("Energy info for: CH3homolysis\n------------------------------------------------------------\n")
		oF.write("Electronic energy (Def2-SVP) = " + ch3_Def2_SVP_lib["homolysis"] + "\n")
		oF.write("Electronic energy (Def2-TZVPD) = " + ch3_Def2_TZVPD_lib["homolysis"] + "\n")
		# oF.write("Zero-point energy = " + ch3_homolysis_lib["zero-point"] + "\n")
		# oF.write("Enthalpy energy = " + ch3_homolysis_lib["enthalpy"] + "\n")
		# oF.write("Free energy = " + ch3_homolysis_lib["free"] + "\n")
	oF.write("------------------------------------------------------------\n")
	if "heterolysisNegative" in keyWords:
		oF.write("Energy info for: CH3(-)heterolysis\n------------------------------------------------------------\n")
		oF.write("Electronic energy (Def2-SVP) = " + ch3_Def2_SVP_lib["heterolysisNegativeMetal"] + "\n")
		oF.write("Electronic energy (Def2-TZVPD) = " + ch3_Def2_TZVPD_lib["heterolysisNegativeMetal"] + "\n")
		# oF.write("Zero-point energy = " + ch3_heterolysisPositiveMetal_lib["zero-point"] + "\n")
		# oF.write("Enthalpy energy = " + ch3_heterolysisPositiveMetal_lib["enthalpy"] + "\n")
		# oF.write("Free energy = " + ch3_heterolysisPositiveMetal_lib["free"] + "\n")
	oF.write("------------------------------------------------------------\n")
	if "heterolysisPositive" in keyWords:
		oF.write("Energy info for: CH3(+)heterolysis\n------------------------------------------------------------\n")
		oF.write("Electronic energy (Def2-SVP) = " + ch3_Def2_SVP_lib["heterolysisPositiveMetal"] + "\n")
		oF.write("Electronic energy (Def2-TZVPD) = " + ch3_Def2_TZVPD_lib["heterolysisPositiveMetal"] + "\n")
		# oF.write("Zero-point energy = " + ch3_heterolysisNegativeMetal_lib["zero-point"] + "\n")
		# oF.write("Enthalpy energy = " + ch3_heterolysisNegativeMetal_lib["enthalpy"] + "\n")
		# oF.write("Free energy = " + ch3_heterolysisNegativeMetal_lib["free"] + "\n")
	oF.write("------------------------------------------------------------\n\n")
	oF.write("Energy info for: " + metal_MetalAndCh3_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
	oF.write("Electronic energy (Def2-TZVPD) = " + metal_MetalAndCh3_lib["ELECTRONIC_ENERGY"] + "\n")
	# oF.write("Zero-point energy = " + metal_MetalAndCh3_lib["ZERO_POINT_ENERGY"] + "\n")
	# oF.write("Enthalpy energy = " + metal_MetalAndCh3_lib["THERMAL_ENTHALPY"] + "\n")
	# oF.write("Free energy = " + metal_MetalAndCh3_lib["FREE_ENERGY"] + "\n")
	oF.write("------------------------------------------------------------\n")
	if "homolysis" in keyWords:
		oF.write("Energy info for: " + metal_homolysis_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
		oF.write("Electronic energy (Def2-TZVPD) = " + metal_homolysis_lib["ELECTRONIC_ENERGY"] + "\n")
		# oF.write("Zero-point energy = " + metal_homolysis_lib["ZERO_POINT_ENERGY"] + "\n")
		# oF.write("Enthalpy energy = " + metal_homolysis_lib["THERMAL_ENTHALPY"] + "\n")
		# oF.write("Free energy = " + metal_homolysis_lib["FREE_ENERGY"] + "\n")
	oF.write("------------------------------------------------------------\n")
	if "heterolysisPositive" in keyWords:
		oF.write("Energy info for: (M+) " + metal_heterolysisPositive_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
		oF.write("Electronic energy (Def2-TZVPD) = " + metal_heterolysisPositive_lib["ELECTRONIC_ENERGY"] + "\n")
		# oF.write("Zero-point energy = " + metal_heterolysisPositive_lib["ZERO_POINT_ENERGY"] + "\n")
		# oF.write("Enthalpy energy = " + metal_heterolysisPositive_lib["THERMAL_ENTHALPY"] + "\n")
		# oF.write("Free energy = " + metal_heterolysisPositive_lib["FREE_ENERGY"] + "\n")
	oF.write("------------------------------------------------------------\n")
	if "heterolysisNegative" in keyWords:
		oF.write("Energy info for: (M-) " + metal_heterolysisNegative_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
		oF.write("Electronic energy (Def2-TZVPD) = " + metal_heterolysisNegative_lib["ELECTRONIC_ENERGY"] + "\n")
		# oF.write("Zero-point energy = " + metal_heterolysisNegative_lib["ZERO_POINT_ENERGY"] + "\n")
		# oF.write("Enthalpy energy = " + metal_heterolysisNegative_lib["THERMAL_ENTHALPY"] + "\n")
		# oF.write("Free energy = " + metal_heterolysisNegative_lib["FREE_ENERGY"] + "\n")
	oF.write("------------------------------------------------------------\n\n")
	if "homolysis" in keyWords:
		if "homolysis" in sucessfulRunTypes:
			oF.write("Bond Energies: (M-Homolysis) " + metal_homolysis_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
			oF.write("Electronic bond energy (Def2-SVP): " + str(metal_homolysis_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY"]) + " = " + str(float(metal_homolysis_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY"]) * 627.5) + " kcal/mol\n")
			oF.write("Electronic bond energy (Def2-TZVPD): " + str(metal_homolysis_lib["ELECTRONIC_BOND_ENERGY"]) + " = " + str(float(metal_homolysis_lib["ELECTRONIC_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Zero-point bond energy: " + str(metal_homolysis_lib["ZERO_POINT_BOND_ENERGY"]) + " = " + str(float(metal_homolysis_lib["ZERO_POINT_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Enthalpy: " + str(metal_homolysis_lib["ENTHALPY_BOND_ENERGY"]) + " = " + str(float(metal_homolysis_lib["ENTHALPY_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Free bond energy: " + str(metal_homolysis_lib["FREE_BOND_ENERGY"]) + " = " + str(float(metal_homolysis_lib["FREE_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
		if "homolysis-singlePoint" in sucessfulRunTypes:
			oF.write("Single-Point electronic bond energy (Def2-SVP): " + str(metal_homolysis_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"]) + " = " + str(float(metal_homolysis_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"]) * 627.5) + " kcal/mol\n")
			oF.write("Single Point electronic bond energy (Def2-TZVPD): " + str(metal_homolysis_lib["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"]) + " = " + str(float(metal_homolysis_lib["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"]) * 627.5) + " kcal/mol\n")
	oF.write("------------------------------------------------------------\n")
	if "heterolysisPositive" in keyWords:
		if "heterolysisPositive" in sucessfulRunTypes:
			oF.write("Bond Energies: ((M+) Heterolysis) " + metal_heterolysisPositive_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
			oF.write("Electronic bond energy (Def2-SVP): " + str(metal_heterolysisPositive_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY"]) + " = " + str(float(metal_heterolysisPositive_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY"]) * 627.5) + " kcal/mol\n")
			oF.write("Electronic bond energy (Def2-TZVPD): " + str(metal_heterolysisPositive_lib["ELECTRONIC_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisPositive_lib["ELECTRONIC_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Zero-point bond energy: " + str(metal_heterolysisPositive_lib["ZERO_POINT_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisPositive_lib["ZERO_POINT_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Enthalpy: " + str(metal_heterolysisPositive_lib["ENTHALPY_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisPositive_lib["ENTHALPY_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Free bond energy: " + str(metal_heterolysisPositive_lib["FREE_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisPositive_lib["FREE_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
		if "heterolysisPositive-singlePoint" in sucessfulRunTypes:
			oF.write("Single-Point electronic bond energy (Def2-SVP): " + str(metal_heterolysisPositive_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"]) + " = " + str(float(metal_heterolysisPositive_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"]) * 627.5) + " kcal/mol\n")
			oF.write("Single Point electronic bond energy (Def2-TZVPD): " + str(metal_heterolysisPositive_lib["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"]) + " = " + str(float(metal_heterolysisPositive_lib["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"]) * 627.5) + " kcal/mol\n")
	oF.write("------------------------------------------------------------\n")
	if "heterolysisNegative" in keyWords:
		if "heterolysisNegative" in sucessfulRunTypes:
			oF.write("Bond Energies: ((M-) Heterolysis) " + metal_heterolysisNegative_lib["FILE_NAME"] + "\n------------------------------------------------------------\n")
			oF.write("Electronic bond energy (Def2-SVP): " + str(metal_heterolysisNegative_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY"]) + " = " + str(float(metal_heterolysisNegative_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY"]) * 627.5) + " kcal/mol\n")
			oF.write("Electronic bond energy (Def2-TZVPD): " + str(metal_heterolysisNegative_lib["ELECTRONIC_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisNegative_lib["ELECTRONIC_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Zero-point bond energy: " + str(metal_heterolysisNegative_lib["ZERO_POINT_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisNegative_lib["ZERO_POINT_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Enthalpy: " + str(metal_heterolysisNegative_lib["ENTHALPY_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisNegative_lib["ENTHALPY_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
			# oF.write("Free bond energy: " + str(metal_heterolysisNegative_lib["FREE_BOND_ENERGY"]) + " = " + str(float(metal_heterolysisNegative_lib["FREE_BOND_ENERGY"]) * 627.5) + " kcal/mol\n")
		if "heterolysisNegative-singlePoint" in sucessfulRunTypes:
			oF.write("Single-Point electronic bond energy (Def2-SVP): " + str(metal_heterolysisNegative_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"]) + " = " + str(float(metal_heterolysisNegative_lib["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"]) * 627.5) + " kcal/mol\n")
			oF.write("Single Point electronic bond energy (Def2-TZVPD): " + str(metal_heterolysisNegative_lib["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"]) + " = " + str(float(metal_heterolysisNegative_lib["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"]) * 627.5) + " kcal/mol\n")
	oF.write("------------------------------------------------------------\n")
	oF.close()
	
#Creates the completed output file
def finalOutput(sucessfulRunTypes):
	if os.path.isdir("../completed"):
		os.chdir("../completed")
	elif os.path.isdir("./completed"):
		os.chdir("./completed")
	elif os.path.isdir("completed/"):
		os.chdir("completed/")
	divider = "-" * 50 + '\n'
	results = []
	conformers_found = 0
	if not os.path.isfile(str(os.environ.get("MOLECULE_NAME")) + ".out"):
		oF = open(str(os.environ.get("MOLECULE_NAME")) + ".out", 'w')
	else:
		oF = open(str(os.environ.get("MOLECULE_NAME")) + ".out", "a")
	oF.write("\nBraden Bond Output File\n")
	time = os.getenv("TIME")
	# oF.write("Trying to write time to final output file\n")
	oF.write("Time taken: " + time + "\n")
	if "True" not in os.getenv("NO_ERRORS"):
		oF.write("Process incomplete because: " + os.getenv("NO_ERRORS"))
		oF.write("\n")
		oF.close()
		return
	# oF.write("After trying to write time to final output file\n")
	oF.write("\n Energy Summaries\n------------------------------------------------------------\n")
	oF.close()
	for file in os.listdir(os.getcwd()):
		results.append(outputFunc(file))
		if file.split(".")[1] == "log":
			if "homolysis" or "heterolysis" in file:
				extractEnergy(file, str(os.environ.get("MOLECULE_NAME")) + ".out")
	writeOutputSummaryFile(str(os.environ.get("MOLECULE_NAME")) + ".out", sucessfulRunTypes)
	# oF.close()

#will convert the xyz format that Billy and Zack's code outputs, to the xyz format that BradenBond expects
def convertFromBillyXyz(xyz_file):
	moleculeName = (xyz_file.split("/")[-1]).split(".")[0]
	os.environ['MOLECULE_NAME'] = str(moleculeName)
	iF = open(xyz_file)
	atom_count = iF.readline()
	line = iF.readline()
	charge = ((line.split("|")[1]).split("= ")[1]).split()[0]
	csd_code = ((line.split("|")[0]).split("= ")[1]).split()[0]
	carbon_index = ((line.split("|")[5]).split(": ")[1]).split()[0]
	second_line_to_write = "F:0-0;P:" + charge + ",1;CI:" + carbon_index +";"+ csd_code + ";\n"
	os.environ["CSD_CODE"] = str(csd_code)
	geometry_data = iF.readlines()
	oF = open(xyz_file, "w")
	oF.write(atom_count)
	oF.write(second_line_to_write)
	for line in geometry_data:
		oF.write(line)
	oF.close()

#first part of runCrest separated in case user runs without crest ??? optimize later
#sets the charge multiplicity and bonds in the input array
#VERSION 2
def setChargeMultBonds(xyz_file, leniency, inputs):
	testFile = open("envVarTestFile", "w")
	bonds = []
	libFile = open(xyz_file, "r")
	header = libFile.readline()
	bonds_line = libFile.readline() #get rid of F: in bond freeze list
	if "|" in bonds_line:
		libFile.close()
		convertFromBillyXyz(xyz_file)
		libFile = open(xyz_file, "r")
		header = libFile.readline()
		bonds_line = libFile.readline()
	else:
		if os.environ.get('CSD_CODE') is None:
			csdCode = bonds_line.split(";")[3]
			moleculeName = xyz_file.split("/")[-1]
			moleculeName = moleculeName.split(".")[0]
			testFile.write("Molecule name as saved: " + str(moleculeName) + "\n")
			os.environ['CSD_CODE'] = str(csdCode)
			os.environ['MOLECULE_NAME'] = str(moleculeName)
	bond_strings = (bonds_line.split(';')[0]).split(':')[1]
	testFile.write("Bond Strings: " + bond_strings + "\n")
	charge_multiplicity = (bonds_line.split(';')[1]).split(":")[1]
	inputs["charge"] = charge_multiplicity.split(',')[0]
	inputs["spin"] = charge_multiplicity.split(',')[1]
	testFile.write("Charge: " + inputs["charge"] + "\n")
	testFile.write("Spin: " + inputs["spin"] + "\n")
	testFile.write("argv[2]: " + str(sys.argv[2]) + "\n")
	testFile.write("argv[1]: " + str(sys.argv[1]) + "\n")
	header += "\n"
	if "CI" in bonds_line: 
		carbon_index = (bonds_line.split(";")[2]).split(":")[1]
		testFile.write("Bonds line info:" + bonds_line)
		testFile.write("Carbon index info: " + carbon_index + "\n")
		os.environ["CARBON_INDEX"] = str(carbon_index)
	if bond_strings[0] == '':
		you_found_the_easter_egg = 5
	else:
		#for bond in bond_strings:
		atoms = bond_strings.split('-')
		bonds.append([str(int(atoms[0]) + 1), str(int(atoms[1]) + 1)]) #chem programs are 1 based so you need to add one
	coords = libFile.readlines()
	inputs["bonds"] = bonds
	libFile.close()

def runCrest(xyz_file, leniency, inputs):
	# print("in runCrest")
	coords = []
	header = ''
	bonds = []
	libFile = open(xyz_file, "r")
	header = libFile.readline()
	bonds_line = libFile.readline()[2:] #get rid of F: in bond freeze list
	bond_strings = bonds_line.split(';')
	bond_strings.pop()
	charge_multiplicity = bond_strings.pop()
	charge_multiplicity = charge_multiplicity[2:]
	charge_multiplicity = charge_multiplicity.split(',')
	inputs["charge"]=charge_multiplicity[0]
	inputs["spin"]=charge_multiplicity[1]
	#with open(xyz_file, "r") as coorid_file:
	header += "\n"
	bond_strings = bond_strings[0].split(',')
	if bond_strings[0] == '':
		you_found_the_easter_egg = 5
	else:
		for bond in bond_strings:
			atoms = bond.split('-')
			bonds.append([str(int(atoms[0]) + 1), str(int(atoms[1]) + 1)]) #chem programs are 1 based so you need to add one
	coords = libFile.readlines()
	inputs["bonds"] = bonds
	os.chdir("crest")
	with open("cinp", "w") as constraint_file:
		constraint_file.write("$constrain\n")
		constraint_file.write("force constant = 0.5\n")
		# constraint_file.write("reference=coords.xyz\n")
		for atoms in bonds:
			constraint_file.write("  distance: " + atoms[0] + ", " + atoms[1] + ", auto\n")
		constraint_file.write("$chrg " + inputs["charge"] + "\n")
		constraint_file.write("$spin " + str(int(inputs["spin"])-1) + "\n") #why -1 on the spin??
		constraint_file.write("$end\n")
	with open("coords.xyz", "w") as crest_coords:
		crest_coords.write(header)
		crest_coords.writelines(coords)
	concatList = []
	if "temperature" in inputs:
		concatList.add("-mdtemp")
		concatList.add(inputs['temperature'])
	if "solvent" in inputs:
		concatList.add("-g")
		concatList.add(inputs['solvent'])
	crest_file = os.getenv('CREST_FILE') #??? where is this env variable set
	run_crest = subprocess.Popen([crest_file, "coords.xyz", "-cinp", "cinp", "--noreftopo", "-ewin", leniency] + concatList) # + [">","crest.out"]) ??? norefttopo flag will
	run_crest.wait()
	try:
		shutil.copy("crest_conformers.xyz", "../crest_conformers.xyz")
		os.chdir("../")
		tempFile = open("crest_conformers.xyz","a")
	except Exception as e:
		print(e)
		sys.exit()
	origFile = open(xyz_file,"r")
	line = origFile.readline()
	#adds the original xyz file into the crest_conformers.xyz, check with someone to see if crest already puts this in???
	#does crest_conformers.xyz order by quality or is it random???
	while line: 
		if "F:" in line:
			tempFile.write("0\n")
		else:
			tempFile.write(line)
		line = origFile.readline()
	tempFile.close()
	origFile.close()

#########################################HOMOLYSIS AND HETEROLYSIS FUNCTIONS AND LIBRARIES###########################################
charge_dict = {
    "homolysis": 0,
    "heterolysisNegativeMetal": -1,
    "heterolysisPositiveMetal": 1
}

mult_dict = {
	"homolysis": 1,
    "heterolysisNegativeMetal": 0,
    "heterolysisPositiveMetal": 0
}

#this function will write the various energies of the ch3 group
def writeCh3Energy(type, final_output_file):
	oF = open(final_output_file, "a")
	if type == "heterolysisNegativeMetal":
		electronicEnergy = ch3_heterolysisNegativeMetal_lib["electronic"]
		zeroPointEnergy = ch3_heterolysisNegativeMetal_lib["zero-point"]
		enthalpy = ch3_heterolysisNegativeMetal_lib["enthalpy"]
		free = ch3_heterolysisNegativeMetal_lib["free"]
	elif type == "heterolysisPositiveMetal":
		electronicEnergy = ch3_heterolysisPositiveMetal_lib["electronic"]
		zeroPointEnergy = ch3_heterolysisPositiveMetal_lib["zero-point"]
		enthalpy = ch3_heterolysisPositiveMetal_lib["enthalpy"]
		free = ch3_heterolysisPositiveMetal_lib["free"]
	elif type == "homolysis":
		electronicEnergy = ch3_homolysis_lib["electronic"]
		zeroPointEnergy = ch3_homolysis_lib["zero-point"]
		enthalpy = ch3_homolysis_lib["enthalpy"]
		free = ch3_homolysis_lib["free"]
	# oF.write("------------------------------------------------------------\n" + "Energy info for ch3 group for " + type + " :\n\n")
	oF.write("Ch3 Electronic energy = " + electronicEnergy + "\n")
	os.environ["CH3_ELECTRONIC"] = electronicEnergy
	oF.write("Ch3 sum of electronic and zero point energies = " + zeroPointEnergy + "\n")
	os.environ["CH3_ZERO_POINT"] = zeroPointEnergy
	oF.write("Ch3 sum of electronic and thermal enthalpies = " + enthalpy + "\n")
	os.environ["CH3_ENTHALPY"] = enthalpy
	oF.write("Ch3 free energy = " + free + "\n")
	os.environ["CH3_FREE"] = free
	# oF.write("------------------------------------------------------------\n")
	oF.close()

#this function will parse the given log file and output various energies and write them to the final output file
def extractEnergy(log_file, final_out_file):
	firstLog = False
	singlePointBool = False
	if "heterolysisNegative" in log_file:
		if "singlePoint" in log_file:
			singlePointBool = True
		type = "heterolysisNegativeMetal"
		dictName = metal_heterolysisNegative_lib
		dictName["FILE_NAME"] = log_file
		# ch3DictName = ch3_heterolysisNegativeMetal_lib
	elif "heterolysisPositive" in log_file:
		if "singlePoint" in log_file:
			singlePointBool = True
		type = "heterolysisPositiveMetal"
		dictName = metal_heterolysisPositive_lib
		dictName["FILE_NAME"] = log_file
		# ch3DictName = ch3_heterolysisPositiveMetal_lib
	elif "homolysis" in log_file:
		if "singlePoint" in log_file:
			singlePointBool = True
		type = "homolysis"
		dictName = metal_homolysis_lib
		dictName["FILE_NAME"] = log_file
		# ch3DictName = ch3_homolysis_lib
	else:
		firstLog = True #this will set the bond energy info for the first log file as evioronment variables
		dictName = metal_MetalAndCh3_lib
		dictName["FILE_NAME"] = log_file
	iF = open(log_file, "r")
	energy = 'not found'
	second_last_energy = 'not found'
	both = False #this specifies wether or not we are using two basis sets, in which case we need to save two energies (the last and second-to-last)
	for line in iF:
		if "SCF Done:" in line:
			second_last_energy = energy
			energy = line.split()[4]
	if singlePointBool:
		dictName["ELECTRONIC_ENERGY_SINGLE_POINT"] = energy
	else:
		dictName["ELECTRONIC_ENERGY"] = energy
	if firstLog:
		os.environ["FIRST_LOG_ELECTRONIC_DEF2SVP"] = str(second_last_energy)
		os.environ["FIRST_LOG_ELECTRONIC_DEF2TZVPD"] = str(energy)
	if not firstLog:
		#bond energy is ch3 energy + metalw/out ch3 energy - original log energy
		#Kcal/mol is the bond energy*627.5
		electronicBondEnergy = (float(ch3_Def2_TZVPD_lib[type]) + float(energy)) - float(os.environ["FIRST_LOG_ELECTRONIC_DEF2TZVPD"])
		if singlePointBool:
			dictName["ELECTRONIC_BOND_ENERGY_SINGLE_POINT"] = electronicBondEnergy
		else:
			dictName["ELECTRONIC_BOND_ENERGY"] = electronicBondEnergy
		secondToLastElectronicBondEnergy = (float(ch3_Def2_SVP_lib[type]) + float(second_last_energy)) - float(os.environ["FIRST_LOG_ELECTRONIC_DEF2SVP"])
		if singlePointBool:
			dictName["SECOND_TO_LAST_ELECTRONIC_ENERGY_SINGLE_POINT"] = secondToLastElectronicBondEnergy
		else:
			dictName["SECOND_TO_LAST_ELECTRONIC_ENERGY"] = secondToLastElectronicBondEnergy
	iF.close()

#This function extracts the geometry from the inital optimized log file to use for the hetrolysis and homolysis xyz files.
def parseLog(log_file, new_xyz_file, numAtomsInt, carbonIndex):
	iF = open(log_file, "r")
	line = iF.readline()
	oF = open(new_xyz_file, "a")
	while "Structure from the checkpoint file:" not in line:#line.split()[0].isnumeric():
		# if "A.U. after    1 cycles" in line:
		# 	extractEnergy(iF)
		line = iF.readline()
		if not line:
			# print("end of file reached, line not found\n")
			return
		continue
	counter = 1
	iF.readline()
	iF.readline()
	# iF.readline()
	# iF.readline()
	# iF.readline()
	line = iF.readline()
	for lines_read in range(numAtomsInt + 1):
		if "," in line:
			lineData = line.split(",")
			if counter != carbonIndex:
				symbol = lineData[0]
				# coords = line.split()[2] + " " + line.split()[3] + " " + line.split()[4]
				coords = lineData[2] + " " + lineData[3] + " " + lineData[4]
				oF.write(symbol + " " + coords) #+ "\n")
				line = iF.readline()
				counter += 1
			else:
				line = iF.readline()
				counter += 1 
				continue
	iF.close()
	oF.close()

#gets the necessary information from the header of the original xyz file and modifies it for the new xyz files, then calls
#the parseLog function to greab the actual geometries that fill the new xyz files
def splitMolecule(xyz_file, type_split):
	iF = open(xyz_file, "r")
	first_line = iF.readline()
	orig_charge_line = iF.readline() 
	original_charge = int(((orig_charge_line.split(";")[1]).split(":")[1]).split(",")[0])
	original_multiplicity = int(((orig_charge_line.split(";")[1]).split(":")[1]).split(",")[1])
	metal_charge = original_charge + charge_dict[type_split]
	# metal_multiplicity = original_multiplicity + mult_dict[type_split]
	if "homolysis" in type_split:
		metal_multiplicity = original_multiplicity + 1
	else :
		metal_multiplicity = 1
	# iF = open(xyz_file, "r")
	oF = open("./newMetalXYZs/" + type_split + "-" + xyz_file, "w")
	numAtomsInt = int(first_line.split()[0]) - 4 #because we will be deleting the carbon and the 3 hydrogens
	# second_line = iF.readline()
	if "CI:" in orig_charge_line:
		carbonIndex = int((orig_charge_line.split(";")[2]).split(":")[1])
	else:
		carbonIndex = -1
	oF.write(str(numAtomsInt) + "\n")
	oF.write("F:0-0;" + "P:" + str(metal_charge) + "," + str(metal_multiplicity) + ";" + str(os.environ.get("CSD_CODE")) + ";\n")
	counter = 1
	# for line in range(numAtomsInt + 1):
	# 	line = iF.readline()
	# 	if counter != carbonIndex:
	# 			oF.write(line)
	# 	counter += 1
	iF.close()
	oF.close()
	xyz_path = "./newMetalXYZs/" + type_split + "-" + xyz_file
	log_path = "./modred/" + "conf" + xyz_file.split(".")[0] + ".log"
	parseLog(log_path, xyz_path, numAtomsInt, carbonIndex)

#main helper function for the homolysis and heterolysis calculations 
def homolysisHeterolysis():
	if not os.path.isdir("./newMetalXYZs"):
		os.mkdir("./newMetalXYZs")
	for filename in os.listdir("./"):
		if filename == sys.argv[2].split("/")[-1]:
			xyz_file = filename
			break
	splitMolecule(xyz_file, "homolysis")
	splitMolecule(xyz_file, "heterolysisNegativeMetal")
	splitMolecule(xyz_file, "heterolysisPositiveMetal")
