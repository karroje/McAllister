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
        self.weightedD = eD
        self.hmmChange = 0
        self.countMatrix = M
        if(M == None):
            self.countMatrix = RptMatrixMod.zeros((4,4))
    
    def calculateChange(self, rate=.6571): # default rate come from regression line in regression tab of DistToProb.xlsm
        #Example .6 calculated and .7 expected means more change
        # this implies a negative number so Hmm should have lower percentages
        return ((self.weightedD/self.numRepeatBases) - self.calculatedD)*rate*-1
    
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
        
    def createPartitions(self, repBasesPerPart = 50000, maxGap = 500000, maxBases = 2000000):
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
                if(i.start - prevEnd > maxGap):
                    if(len(self.partitionList) > 0 and 
                       self.currentPartition.startIndex - self.partitionList[-1].endIndex < maxGap):
                        self.partitionList[-1] = self.combinePartitions(self.partitionList[-1], self.currentPartition)
                    # reset partition
                    self.currentPartition = Partition()
                
                self.currentPartition.endIndex = i.finish
                self.currentPartition.countMatrix +=i.M
                self.currentPartition.numRepeatBases += i.finish - i.start
                # Add to the numerator of what will be a weighted value
                self.currentPartition.weightedD += (i.finish - i.start)*self.familiesDict[i.class_name]
                prevEnd = i.finish
                # If number of bases is achieved create the partition
                if(self.currentPartition.numRepeatBases > repBasesPerPart):
                    self.finalizePartition(self.currentPartition)
                    self.partitionList.append(self.currentPartition)
                    self.currentPartition = Partition()
                    self.currentPartition.startIndex= self.partitionList[-1].endIndex
                    
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
                    
                
                
    # Calculates the distances based on the counts    
    def finalizePartition(self,part):
        part.calculatedD =  compute_d(compute_P(part.countMatrix))
        part.hmmChange = part.calculateChange()
            
        
     
            
        
    
#p = PartitionMaker()
#p.createPartitions()  
    
    