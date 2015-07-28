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
from rate_tools import HMMEditor

class HillClimber:
    def __init__(self, hmmName, chrName, makeDirs =True):
        self.folder = "./HillClimber_" + hmmName.split()[0] + "_" + chrName.split()[0] + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "/"
        if(makeDirs):
            os.makedirs(self.folder)
        self.hmmFolder = self.folder + "Hmms/"
        if(makeDirs):
            os.makedirs(self.hmmFolder)
        self.resultsFolder = self.folder + "Results/"
        if(makeDirs):
            os.makedirs(self.resultsFolder)
        self.hmmBaseName = hmmName
        self.hmmBest = self.copyOriginalHmmFile()
        self.hmmBestAmount = 0
        self.hmmMostBases = 0 
        self.annealingRate = .8
        self.minChange = 1
        self.startChangePercent = 5
        self.seqFile = "./FastaFiles/" + chrName
    
    def copyOriginalHmmFile(self):
        newHmmName = self.makeHMMName(0)
        filePath = "./AllHMMs/" + self.hmmBaseName +".hmm"
        shutil.copy(filePath, self.hmmFolder + newHmmName)
        return newHmmName
    
    # run the HMM through on the 
    def runHMM(self,hmmFileName):  
        hmmFile = self.hmmFolder + hmmFileName
        chrName = self.seqFile.split("/")[-1].split(".")[0]
        hmmName = hmmFileName.split(".")[0]
        resultFile = self.resultsFolder + hmmName + ".res"
        hmmerInterOrig = HmmerInterface.HmmerInfo(self.seqFile, hmmFile, resultFile)
        proc = hmmerInterOrig.findResults()
        proc.wait()
        numBasesFound = self.analyzeResultFile(resultFile)
        return numBasesFound
    # determine number of bases found by the 
    def analyzeResultFile(self, resultFile):
        numBases = 0
        results = FileMaker.parseResultsFile(resultFile, resultFile + "cleaned", 0)
        return results[0]
    
    def updateBestFile(self, hmmName, numFoundBases, amt):
        print("Compare: " + hmmName +" : " + str(numFoundBases) +
              " ; " + self.hmmBest +" : " + str(self.hmmMostBases))
        if(numFoundBases > self.hmmMostBases):
            self.hmmBest = hmmName
            self.hmmMostBases = numFoundBases
            self.hmmBestAmount = amt
        
    def checkIfHMMFileExists(self,hmmFile):
        fileNames = os.listdir(self.hmmFolder)
        return (hmmFile in fileNames)
    def createHMMFile(self,prevFile, adjustment, newHmmName):
        if( not self.checkIfHMMFileExists(newHmmName)):
            prevProfile = HMMEditor.hmmprofile()
            prevProfile.parseFile(self.hmmFolder + prevFile)
            prevProfile.addToDecimalEmissionProbs(adjustment, .95, .01, True)
            prevProfile.convertDecimalProbsToEmission()
            prevProfile.writeFile(self.hmmFolder + newHmmName)
    
    def makeHMMName(self, percentChange):
        return self.hmmBaseName +"_" + str(percentChange) +".hmm"
    
    def createRunUpdate(self, hmmName, diffBestChange, totChgAmt):
        self.createHMMFile(self.hmmBest,(diffBestChange + 0.0)/100,hmmName)
        basesFound = self.runHMM(hmmName)
        self.updateBestFile(hmmName, basesFound, totChgAmt)
    # Main method that runs the trials
    def runHillClimber(self):
        self.createRunUpdate(self.hmmBest, 0, 0)
        change = self.startChangePercent
        while change > self.minChange:
            intChange = int(round(change,0))
            hmmHigherAmt = self.hmmBestAmount + intChange
            hmmLowerAmt = self.hmmBestAmount - intChange
            hmmHigherName = self.makeHMMName(hmmHigherAmt)
            hmmLowerName = self.makeHMMName(hmmLowerAmt)
            if(not self.checkIfHMMFileExists(hmmHigherName)):
                self.createRunUpdate(hmmHigherName, intChange, hmmHigherAmt)
            if(not self.checkIfHMMFileExists(hmmLowerName)):
                self.createRunUpdate(hmmLowerName, -intChange, hmmLowerAmt)
            change *= self.annealingRate
        self.finalCompare()
        print(self.hmmBest)
        print(self.hmmBestAmount)
    def finalCompare(self):
        origList = self.getRepeatList(self.getCleanResultFileName(0))
        bestList = self.getRepeatList(self.getCleanResultFileName(self.hmmBestAmount))
        origList.sort(key=lambda x: x.start)
        bestList.sort(key=lambda x: x.start)
        uncovered = self.countListMisses(origList, bestList)
        uncoveredInOrig = uncovered[0]
        uncoveredInBest = uncovered[1]
        print(uncoveredInOrig)
        print(uncoveredInBest)
        print(self.internalOverlaps(origList))
        print(self.internalOverlaps(bestList))
        
    def countListMisses(self, origList, bestList):
        origUncovered = 0
        bestUncovered = 0 
        origIndex = 0
        bestIndex = 0
        while origIndex < len(origList) and bestIndex < len(bestList):
            if(origList[origIndex].start < bestList[bestIndex].start):
                #Best has to not started by time orig finishes
                if origList[origIndex].end < bestList[bestIndex].start:
                    if not origList[origIndex].covered:
                        origUncovered +=1
                        print ("Orig uncovered " + str(origList[origIndex]))
                    origIndex +=1
                # Best has to have started, but not finished
                elif origList[origIndex].end <= bestList[bestIndex].end:
                    origList[origIndex].covered=True
                    bestList[bestIndex].covered=True
                    origIndex += 1
                # orig end > best end && orig start < best start
                else:
                    origList[origIndex].covered=True
                    bestList[bestIndex].covered=True
                    bestIndex += 1
            else:
                # best start < orig start && best end < orig start
                if origList[origIndex].start > bestList[bestIndex].end:
                    if not bestList[bestIndex].covered:
                        bestUncovered +=1
                        print ("best uncovered " + str(bestList[bestIndex]))
                    bestIndex +=1
                #best start < orig start && best end < orig end
                elif origList[origIndex].end >= bestList[bestIndex].end:
                    origList[origIndex].covered=True
                    bestList[bestIndex].covered=True
                    bestIndex +=1
                else: 
                    origList[origIndex].covered=True
                    bestList[bestIndex].covered=True
                    origIndex +=1
        if origIndex < len(origList):
            for i in range(origIndex, len(origList)):
                if(not origList[i].covered):
                    origUncovered +=1
        if bestIndex < len(bestList):
            for i in range(bestIndex, len(bestList)):
                if(not bestList[i].covered):
                    bestUncovered +=1
            
        return [origUncovered, bestUncovered]
                    
    
    def getRepeatList(self, resultFileName):
        resFileHandle = open(resultFileName,'r')
        repeatList = []
        index = 0
        #Read first 3 header Lines
        for line in resFileHandle:
            lineList = line.split()
            if index > 2:
                repeatList.append(RepeatInstance(int(lineList[4]),int(lineList[5]),lineList[0]))
            index+=1
        return repeatList
    
    def getCleanResultFileName(self, change):
        resultFile = self.resultsFolder + self.hmmBaseName + "_" + str(change) + ".rescleaned"
        return resultFile
    def internalOverlaps(self, list):
        overlaps = 0
        prevEnd = -1000000000
        for repeat in list:
            if(repeat.start  < prevEnd):
                overlaps +=1
            prevEnd = repeat.end
        return overlaps
class RepeatInstance:
    def __init__(self,s, e, sc):
        self.start = s
        self.end = e
        self.score =sc
        self.covered = False
    def __repr__(self):
        return ' {} {} {}'.format(self.start,
                                  self.end,self.score)
    def __str__(self):
        return ' {} {} {}'.format(self.start,
                                  self.end,self.score)
 
    def __cmp__(self, other):
        if hasattr(other, 'start'):
            return self.start.__cmp__(other.start)


if __name__ == "__main__":
    hc = HillClimber("L2", "chr22.fa")
    hc.runHillClimber()
    
    #hc = HillClimber("MLTA", "chr22.fa")  
    #hc.resultsFolder = "./HillClimber_AluSx_chr22.fa2015-07-28 15:27:23/Results/"
    #hc.hmmBestAmount= -3
    #hc.finalCompare()
        