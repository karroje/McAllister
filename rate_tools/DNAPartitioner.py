"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from numpy import rate
import operator
import RptMatrixMod
import operator
import pickle
from asymm_tools import compute_d
from asymm_tools import compute_P

class Partition:
    def __init__(self, s=0, e =0, b=0, cD=0,eD=0, M=None):         
        self.startIndex = s  # Starting index in chromosone
        self.endIndex = e # Ending index in Chromosone
        self.numRepeatBases = b #Number of repeat bases in Partition
        self.calculatedD =  cD #Distance based on counts
        self.calculatedD2 = cD
        self.weightedD = eD
        self.hmmChange = 0
        self.countMatrix = M
        if(M == None):
            self.countMatrix = RptMatrixMod.zeros((4,4))
    
    def calculateChange(self, rate=.6571): # default rate come from regression line in regression tab of DistToProb.xlsm
        #Example .6 calculated and .7 expected means more change
        # this implies a negative number so Hmm should have lower percentages
        return ((self.weightedD/self.numRepeatBases) - self.calculatedD)*rate
        #return ((self.weightedD/self.numRepeatBases) - self.calculatedD2/self.numRepeatBases)*rate*-1
    
class PartitionMaker:
    
    def __init__(self, search= "MIR" ,psm = "/Users/MikeMcAllister/cache/human/hg18/seq/rmsk/chr11FromFileMaker.psm", famFile = "../rpt_list.txt", famDist = "../GlobalIDs.txt"):
        self.searchRepeat = search
        self.psmFile = psm
        self.famFile = famFile
        self.famDistFile = famDist
        self.familiesDict = self.createDict(self.famDistFile)
        self.includedFamDict ={}
        #self.familiesList = self.makeList() 
        self.pruneFamilies(.1)
        self.currentPartition = Partition()
        self.partitionList =[]
    
    def createDict(self, distFile):
        fh = open(distFile, 'r')
        familyDict = {}
        for line in fh:
            line = line.strip()
            list = line.split(" ")
            familyDict[list[0]] = float(list[1])
        return familyDict
    
    #def makeList(self):
    #    return sorted(self.familiesDict.items(), key=operator.itemgetter(1))
   
    def pruneFamilies(self, buffer = 1):
        if (self.searchRepeat in self.familiesDict):
            includedFamDict = {}
            center = self.familiesDict[self.searchRepeat]
            for k,v in self.familiesDict.iteritems():
                if(v > center - buffer and v< center + buffer):
                    self.includedFamDict[k] = v
        else:
            print("Trying to prune with unknown family")
        
    def createPartitions(self, repBasesPerPart = 50000, maxGap = 500000, maxBases = 1000000):
        rpt_list = RptMatrixMod.load_psm(self.psmFile)
        prevEnd = rpt_list[0].start
        self.currentPartition.startIndex = rpt_list[0].start
        for i in rpt_list:
            if i.class_name in self.includedFamDict:
                # Place start index in partition if partition is empty
                if(self.currentPartition.startIndex == 0):
                    self.currentPartition.startIndex = i.start 
                # If the gaps between repeats is too high due to whatever start a new partition and combine the current
                # one with the previous partition
                self.testGap(prevEnd, i, maxGap)
                self.testMaxDistanceBases(i, maxBases)
                self.currentPartition.endIndex = i.finish
                self.currentPartition.countMatrix +=i.M
                self.currentPartition.numRepeatBases += i.finish - i.start + 1
                i.M = self.testNegativeDiagonal(i.M)
                self.currentPartition.calculatedD2 += (i.finish - i.start)*compute_d(compute_P(i.M))
                # Add to the numerator of what will be a weighted value
                self.currentPartition.weightedD += (i.finish - i.start)*self.familiesDict[i.class_name]
                prevEnd = i.finish
                # If number of bases is achieved create the partition
                self.testTotalBases(repBasesPerPart)
        self.combinePartitions(self.partitionList[-1], self.currentPartition,)
        
    def testMaxDistanceBases(self, repeat, maxBases):
        if(repeat.finish - self.currentPartition.startIndex > maxBases):
            self.createNewPartition()
            self.currentPartition.startIndex= repeat.start
            
    
    def endPartition(self):
        #print(self.currentPartition.countMatrix)
        self.finalizePartition(self.currentPartition)
        self.partitionList.append(self.currentPartition)
    def createNewPartition(self):           
        self.currentPartition = Partition()
        self.currentPartition.startIndex= self.partitionList[-1].endIndex + 1
        
    def testTotalBases(self, repBasesPerPart):
        if(self.currentPartition.numRepeatBases > repBasesPerPart):
            self.endPartition()
            self.createNewPartition()
            
    def testGap(self, prevEnd, repeat, maxGap):
        if(repeat.start - prevEnd > maxGap):
            if(len(self.partitionList) > 0 and 
               self.currentPartition.startIndex - self.partitionList[-1].endIndex < maxGap):
                self.partitionList[-1] = self.combinePartitions(self.partitionList[-1], self.currentPartition)
                    # reset partition
                self.createNewPartition()
                    
    def printPartitionList(self):
        index = 1
        for j in self.partitionList:
            print("Partion " + str(index) + ": Start " + str(j.startIndex) + "; End: " + str(j.endIndex)
                  + "; Length: " + str(j.numRepeatBases))
            print("ExpectedD " + str(j.weightedD/j.numRepeatBases) + " ; CalculatedD: " + str(j.calculatedD)
                  + "; Difference: " + str(j.weightedD/j.numRepeatBases-j.calculatedD))
            index +=1
                    
    def combinePartitions(self,p1,p2):
        p1.startIndex = min(p1.startIndex, p2.startIndex)  
        p1.endIndex = max(p1.endIndex, p2.endIndex) 
        p1.weightedD = p1.weightedD + p2.weightedD
        p1.numRepeatBases += p2.numRepeatBases
        p1.countMatrix += p2.countMatrix
        self.finalizePartition(p1)
        
        return p1
                    
    def testNegativeDiagonal(self, matrix):
        if(matrix[0][0] == 0):
           matrix[0][0] = 1 
        if(matrix[1][1] == 0):
           matrix[1][1] = 1    
        if(matrix[2][2] == 0):
           matrix[2][2] = 1 
        if(matrix[3][3] == 0):
           matrix[3][3] = 1
        return matrix       
                
    # Calculates the distances based on the counts    
    def finalizePartition(self,part):
        part.calculatedD =  compute_d(compute_P(part.countMatrix))
        part.hmmChange = part.calculateChange()
            
        
     
            
        
    
#p = PartitionMaker()
#p.createPartitions()  
    
    