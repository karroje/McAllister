"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import HmmerInterface
import FileMaker
import argparse
import datetime
import os.path
import SimulationResultSummary

class FileInfo:
    def __init__(self):
        self.partitions = []
        self.folder = ""
        self.hmmFolder =""
        self.seqFolder =""
        self.hmmStub = ""
        self.fastaStub = ""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Edit to be simulated results
    parser.add_argument('-f', action='store', dest='fastaFile',
                    default = "./FastaFiles/simulation.fa", help='Set Fasta File')
    parser.add_argument('-m', action='store', dest='hmmFile',
                    default = "./AllHMMs/AluSx.hmm", help='Set Orig HMM file')
    parser.add_argument('-p', action='store', dest='HMMERPATH',
                    default = "/usr/local/bin/", help='Path to HMM software')
    res = "./HMMResults" + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') +"/"
    parser.add_argument('-r', action='store', dest='resultsFolder',
                    default = res, help='Where results will be located')
    parser.add_argument('-s', action='store', dest='psmFile',
                    default = "./PSMs/simulation.psm", help='Where PSM file is located')
    results = parser.parse_args()
    
    fInfo = FileInfo()
    
    #Create the Files to run
    FileMaker.createFiles(results.fastaFile, results.hmmFile,
                          results.psmFile, results.resultsFolder, fInfo)
    
    procs = []
    
    # Fork the processes to run 
    for index in xrange(len(fInfo.partitions)):
        part = fInfo.partitions[index]
        hmmFile = fInfo.hmmFolder + fInfo.hmmStub + str(index) 
        seqFile = fInfo.seqFolder + fInfo.fastaStub + str(index)
        resFile = results.resultsFolder + str(index) + ".res"
        hmmerInter = HmmerInterface.HmmerInfo(seqFile, hmmFile, resFile)
        procs.append(hmmerInter.findResults())
        #Run Original Hmm file
        hmmerInterOrig = HmmerInterface.HmmerInfo(seqFile, results.hmmFile, resFile + "orig")
        procs.append(hmmerInterOrig.findResults())
        procs[-2].wait()
        procs[-1].wait()
    
    #Gather results
    repBasesMod = 0
    repBasesOrig = 0
    changeList  = []
    for index in xrange(len(fInfo.partitions)):
        p = procs[index*2]
        part = fInfo.partitions[index]
        resFile = results.resultsFolder + str(index) + ".res"
        #p.wait()
        modResults = FileMaker.parseResultsFile(resFile, resFile + "cleaned", part.startIndex)
        pOrig = procs[index*2 + 1]
        resFileOrig = results.resultsFolder + str(index) + ".resorig"
        #pOrig.wait()
        origResults =FileMaker.parseResultsFile(resFileOrig, resFileOrig + "cleaned", part.startIndex)
        repBasesMod += modResults[0]
        repBasesOrig += origResults[0]
        if(float(modResults[0])/modResults[1] > float(origResults[0])/origResults[1]):
            print( "More results in modified.  Partition: " + str(index))
        elif(float(modResults[0])/modResults[1] < float(origResults[0])/origResults[1]):
            print("Less results in modified.  Partition: " + str(index))
        else:
            print("Same results.  Partition: " + str(index))
        print(str(modResults[0]) + " : " + str(origResults[0]) + " : "  + str(part.calculateChange()))
        changeList.append(part.calculateChange())
    print(str(repBasesMod) + " : " + str(repBasesOrig))
    FileMaker.aggregrateResultFiles(results.resultsFolder, ".rescleaned")   
    FileMaker.aggregrateResultFiles(results.resultsFolder, ".resorigcleaned") 
    SimulationResultSummary.SimulationResultSummary(results.resultsFolder + "Cumulative.rescleaned", 
                                                    results.resultsFolder + "Cumulative.resorigcleaned",
                                                    changeList = changeList) 
    