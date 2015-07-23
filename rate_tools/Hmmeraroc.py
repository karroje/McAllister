"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import shutil
import HmmerInterface
import FileMaker
import argparse
import datetime
import os
import SimulationResultSummary
import re

class FileInfo:
    def __init__(self):
        self.partitions = []
        self.folder = ""
        self.hmmFolder =""
        self.seqFolder =""
        self.hmmStub = ""
        self.fastaStub = ""


    
# Fork the processes to run 
def runPartitions(fInfo, results, procs):
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
    
def deleteTempFileFolder(fInfo):
    shutil.rmtree(fInfo.folder)
    
def deleteInitialResultFiles(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))

    #Gather results
def summarizeResults(fInfo, results, procs):
    partResults = open(results.resultsFolder + "basesByPartition.txt",'w')
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
            partResults.write( "\nMore results in modified.  Partition: " + str(index))
        elif(float(modResults[0])/modResults[1] < float(origResults[0])/origResults[1]):
            partResults.write("\nLess results in modified.  Partition: " + str(index))
        else:
            partResults.write("\nSame results.  Partition: " + str(index))
        partResults.write("\n" + str(modResults[0]) + " : " + str(origResults[0]) + " : "  + str(part.calculateChange()))
        changeList.append(part.calculateChange())
    partResults.write("\n" + str(repBasesMod) + " : " + str(repBasesOrig))
    FileMaker.aggregrateResultFiles(results.resultsFolder, ".rescleaned")   
    FileMaker.aggregrateResultFiles(results.resultsFolder, ".resorigcleaned") 
    if(results.s):
        SimulationResultSummary.SimulationResultSummary(results.resultsFolder + "Cumulative.rescleaned", 
                                                    results.resultsFolder + "Cumulative.resorigcleaned",
                                                    changeList = changeList)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Edit to be simulated results
    parser.add_argument('-f', action='store', dest='fastaFile',
                    default = "./FastaFiles/simulation.fa", help='Set Fasta File')
    parser.add_argument('-m', action='store', dest='hmmFile',
                    default = "./AllHMMs/MIR.hmm", help='Set Orig HMM file')
    parser.add_argument('-b', action='store', dest='HMMERPATH',
                    default = "/usr/local/bin/", help='Path to HMM software')
    res = "./HMMResults" + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') +"/"
    parser.add_argument('-r', action='store', dest='resultsFolder',
                    default = res, help='Where results will be located')
    parser.add_argument('-p', action='store', dest='psmFile',
                    default = "./PSMs/simulation.psm", help='Where PSM file is located')
    parser.add_argument('-s', action='store_true', help='Run is Simulation.  Make Simulation ' +
                        ' summary files')
    parser.add_argument('-n', action='store_true', help='Make only negative changes to  ' +
                        'original base probabilities')
    parser.add_argument('-d', action='store_false', help='deletes folder used to hold HMMs ' + 
                        ' and fasta files for each partition')
    parser.add_argument('-e', action='store_false', help='deletes original results files ' + 
                        ' and only keeps the cumulative ones')
    results = parser.parse_args()
    print("start: " + str(datetime.datetime.now().time()))
    fInfo = FileInfo()
    #Create the Files to run
    FileMaker.createFiles(results.fastaFile, results.hmmFile,
                          results.psmFile, results.resultsFolder, fInfo, results.n)
    print("Start HMMER: " +str(datetime.datetime.now().time()))
    procs = []
    runPartitions(fInfo, results, procs)
    summarizeResults(fInfo, results, procs)  
    if(results.d):
        deleteTempFileFolder(fInfo)
    if(results.e):
        deleteInitialResultFiles(results.resultsFolder,"^\d+\.(.)*")
    print("End: " + str(datetime.datetime.now().time()))