"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
class partitionStats:
    def __init__(self):
        self.start = 0
        self.end = 0
        self.index = 0
        self.repeatBasesPresent = 0
        self.repeatBasesMatched = 0
        self.repeatBasesIdentified = 0

class SimulationResultSummary:
    def __init__(self, resultModFile, resultOrigFile,
                 rlf = "simulationRepeatIndices.csv",
                 changeList = None):
        self.repeatLocationsFile = rlf
        self.resModFile = resultModFile
        self.resOrigFile = resultOrigFile
        self.repeatLocationList = self.populateRepeatList()
        self.filePath = os.path.split(self.resModFile)[0]
        self.summaryModFile = self.filePath + "/summaryOfMod.txt"
        self.summaryOrigFile = self.filePath + "/summaryOfOrig.txt"
        self.summaryFile = self.filePath + "/summary.txt"
        self.partitionChangeList = changeList
        self.parseResultsFile(self.resModFile, self.summaryModFile)
        self.parseResultsFile(self.resOrigFile, self.summaryOrigFile)
        self.makeCombinedFile()
    
    def populateRepeatList(self):
        repeatFile = open(self.repeatLocationsFile,"r")
        repeatList = []
        for line in repeatFile:
            linePieces = line.split(",")
            start = int(linePieces[0])
            end = int(linePieces[1])
            repeatTuple = (start, end)
            repeatList.append(repeatTuple)
        return repeatList
    
    """
    looks for closest lower starting repeat list index to supplied starting coord
    Exception if repeat start coord is within 5 bases of supplied coord then in 
    does not have to be lower
    """
    def findClosestIndex(self, startCoord):
        beginSearch = 0
        endSearch = len(self.repeatLocationList) - 1
        lastLowerVal = 0
        lastLowerIndex = 0
        unmatched = True
        while(beginSearch <= endSearch and unmatched):
            midIndex = int((beginSearch + endSearch) / 2)
            midIndexStart = self.repeatLocationList[midIndex][0]
            if(abs(startCoord - midIndexStart) < 400):
                return midIndex
            elif(startCoord > midIndexStart):
                lastLowerVal = midIndexStart
                lastLowerIndex = midIndex
                beginSearch = midIndex + 1
            else:
                endSearch = midIndex - 1
        return lastLowerIndex
    
    def overLapCount(self, startCoord, endCoord, repeatIndex):
        repeat = self.repeatLocationList[repeatIndex]
        repeatStart = repeat[0]
        repeatEnd = repeat[1]
        matchedBases = 0
        if (repeatEnd < startCoord):
            return matchedBases
        elif(endCoord < repeatStart):
            print("This line should not be hit check findClosestIndexMethod")
            return matchedBases
        else:
            highStart = max(startCoord,repeatStart)
            lowEnd = min(endCoord, repeatEnd)
            return lowEnd - highStart
    """
    self.start = 0
        self.end = 0
        self.repeatBasesPresent = 0
        self.repeatBasesMatched = 0
        self.repeatBasesIdentified = 0
    """
        
    def writePartitionSummary(self, summaryHandle, partStats):
        summaryHandle.write("Partition " + str(partStats.index)
                            + ": Start " + str(partStats.start) +
                            "; End: " + str(partStats.end))
        
        if(self.partitionChangeList != None):
            summaryHandle.write(" Change to HMM Probs: " 
                                + str(self.partitionChangeList[partStats.index])
                                + "\n")
        
        summaryHandle.write("\tBases Matched: " + str(partStats.repeatBasesMatched) + "\n")
        summaryHandle.write("\tTotal Bases: " + str(partStats.repeatBasesPresent) + "\n")
        summaryHandle.write("\tBases Missed: " + str(partStats.repeatBasesPresent -
                            partStats.repeatBasesMatched) + "\n")
        summaryHandle.write("\tBases Incorrectly identified: " +
                            str(partStats.repeatBasesIdentified - 
                                partStats.repeatBasesMatched)+ "\n")
        
        
    def parseResultsFile(self, resultsFile, summaryFile):
        summaryHandle = open(summaryFile,"w")
        resultsHandle = open(resultsFile,"r")
        prevPartitionStart = -1
        currentPart = partitionStats()
        partIndex = 0
        for line in resultsHandle:
            linePieces = line.split()
            if(prevPartitionStart != int(linePieces[-3])):
                self.writePartitionSummary(summaryHandle, currentPart)
                currentPart = partitionStats()
                currentPart.end = int(linePieces[-1])
                currentPart.start = int(linePieces[-3])
                currentPart.index = partIndex
                prevPartitionStart =  int(linePieces[-3])
                currentPart.repeatBasesPresent = self.countRepeatBasesPresent(currentPart.start, currentPart.end)
                partIndex +=1
            
            start = int(linePieces[4])
            end = int(linePieces[5])
            if(start < end):
                currentPart.repeatBasesIdentified += (end - start) + 1 
                currentPart.repeatBasesMatched += (self.overLapCount(start, end, self.findClosestIndex(start)))
            
    def countRepeatBasesPresent(self, startPartition, endPartition):
        repeatBases  = 0
        for repeat in self.repeatLocationList:
            
            if(repeat[0] > startPartition and repeat[1] < endPartition):
                repeatBases += repeat[1] - repeat[0] + 1
            elif(repeat[0] > startPartition and repeat[1] > endPartition and repeat[0] < endPartition):
                repeatBases += endPartition - repeat[0] + 1
            elif(repeat[0] < startPartition and repeat[1] < endPartition and repeat[1] > startPartition):
                repeatBases += repeat[1] - startPartition + 1
            elif(repeat[0]< startPartition and repeat[1]>endPartition):
                repeatBases += endPartition - startPartition + 1
                
        return repeatBases
    
    def makeCombinedFile(self):
        modFile = open(self.summaryModFile,"r")
        origFile = open(self.summaryOrigFile,"r")
        combinedFile = open(self.summaryFile,"w")
        for modLine in modFile:
            modLine = modLine.strip("\n")
            origLine = origFile.readline().strip("\n")
            if modLine.startswith("Partition"):
                combinedFile.write(modLine)
            else: 
                combinedFile.write(modLine.ljust(40))
                combinedFile.write(" : ")
                combinedFile.write(origLine.ljust(40))
            combinedFile.write("\n")
                
if __name__ == "__main__":
    srs = SimulationResultSummary("HMMResults2015-07-07 16:49:17/Cumulative.rescleaned",
                                  "HMMResults2015-07-07 16:49:17/Cumulative.resorigcleaned")
    
        
        
        
        
                
        