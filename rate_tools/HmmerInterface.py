"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import subprocess
import os



class HmmerInfo:
	def __init__(self, seqFile ="testDNASeq", hmmFile="set", results = "results",
				HMMERPath = "/usr/local/bin/"):
		self.repeat = "ACGTTGCAGCTAGCT"
		self.hmmbuild_location = HMMERPath
		self.repeatAlignFile = "testAlign"
		self.seqFile= seqFile
		if hmmFile == "set":
			self.hmmFile = os.path.join(os.path.dirname(__file__), 'HMMs/testHmm')
		else:
			self.hmmFile = hmmFile
		self.resultRepeatsFile = results
	# needed to chmod + x the binary
	def buildHmm(self):
		#commandList = ["/Users/MikeMcAllister/Google Drive/Grad/Research/hmmer-3.1b1-macosx-intel/binaries/hmmbuild", self.hmmFile, self.repeatAlignFile]
		commandList = [self.hmmbuild_location + "hmmbuild", self.hmmFile, self.repeatAlignFile]
		#subprocess.call("/Users/MikeMcAllister/Google\ Drive/Grad/Research/hmmer-3.1b1-macosx-intel/binaries/hmmbuild --dna " + self.hmmFile + " " + self.repeatAlignFile, shell =True)
		file = open(self.hmmFile,"w")
		file.close()
		subprocess.Popen(commandList)
	def findResults(self):
		#commandList = ["/Users/MikeMcAllister/Google Drive/Grad/Research/hmmer-3.1b1-macosx-intel/binaries/nhmmer", "--dna", self.hmmFile, self.seqFile]
		commandList = [self.hmmbuild_location+"nhmmer", "--dna", self.hmmFile, self.seqFile]
		outFile = open(self.resultRepeatsFile,"w")
		return subprocess.Popen(commandList, stdout = outFile)


if __name__ == "__main__":		
	testHI = HmmerInfo()
	
	#testHI.repeatAlignFile = "testAlignChanges"
	#testHI.buildHmm()
	testHI.findResults()
	testHI.hmmFile = "./HMMs/HMMGreaterChange"
	testHI.resultRepeatsFile ="results2"
	testHI.findResults()