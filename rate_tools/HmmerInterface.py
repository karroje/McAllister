import subprocess
import os

hmmbuild_location = "/usr/local/bin/"

class HmmerInfo:
	def __init__(self):
		self.repeat = "ACGTTGCAGCTAGCT"
		self.repeatAlignFile = "testAlign"
		self.seqFile="testDNASeq"
		self.hmmFile = os.path.join(os.path.dirname(__file__), 'HMMs/testHmm')
		self.resultRepeatsFile = "results"
	# needed to chmod + x the binary
	def buildHmm(self):
		#commandList = ["/Users/MikeMcAllister/Google Drive/Grad/Research/hmmer-3.1b1-macosx-intel/binaries/hmmbuild", self.hmmFile, self.repeatAlignFile]
		commandList = [hmmbuild_location + "hmmbuild", self.hmmFile, self.repeatAlignFile]
		#subprocess.call("/Users/MikeMcAllister/Google\ Drive/Grad/Research/hmmer-3.1b1-macosx-intel/binaries/hmmbuild --dna " + self.hmmFile + " " + self.repeatAlignFile, shell =True)
		file = open(self.hmmFile,"w")
		file.close()
		subprocess.Popen(commandList)
	def findResults(self):
		#commandList = ["/Users/MikeMcAllister/Google Drive/Grad/Research/hmmer-3.1b1-macosx-intel/binaries/nhmmer", "--dna", self.hmmFile, self.seqFile]
		commandList = [hmmbuild_location+"nhmmer", "--dna", self.hmmFile, self.seqFile]
		outFile = open(self.resultRepeatsFile,"w")
		subprocess.Popen(commandList, stdout = outFile)
		
testHI = HmmerInfo()

#testHI.repeatAlignFile = "testAlignChanges"
#testHI.buildHmm()
testHI.findResults()
testHI.hmmFile = "./HMMs/HMMGreaterChange"
testHI.resultRepeatsFile ="results2"
testHI.findResults()