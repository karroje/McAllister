import HMMEditor
import pickle
import RptMatrixMod
import random
import numpy
from Bio.Nexus.Trees import consensus
class SimulationContainer:
    def __init__(self, start = .25, end = .4,incr = .005):
        self.otherRepeatFileNames = ["./HMMs/LTR16a.hmm","./HMMs/MARNA.hmm","./HMMs/MER115.hmm","./HMMs/mlt1l.hmm", "./HMMs/tigger8.hmm"]
        self.repeatHMMs = self.retrieveRepeatList(self.otherRepeatFileNames)
        self.startPercentChange = start
        self.endPercentChange = end
        self.percentIncrement = incr
        self.repeatToFind = self.retrieveHMM("./HMMs/MADE2.hmm")
        self.repeatIndicesFile = "simulationRepeatIndices.csv"
        
        
    def retrieveHMM(self, fileName):
        profile = HMMEditor.hmmprofile()
        profile.parseFile(fileName)
        return profile
    
    def retrieveRepeatList(self,fileNames):
        hmmList = []
        for name in fileNames:
            hmmList.append(self.retrieveHMM(name))
        return hmmList
    
    # Modifies consensus sequence of HMM with probability specified
    def modifySequence(self, hmm, probability):
        modSeq = ""
        for x in range(1,len(hmm.emissionDecimalProbs)):
            consensusIndex = hmm.getConsensusIndex(x)
            origLetter = hmm.alphabet[consensusIndex]
            nextLetter = origLetter
            if random.random() < probability:
                probs = hmm.emissionDecimalProbs[x]
                sumProbs = 0
                
                for i in range(len(probs)):
                    if (i != consensusIndex):
                        sumProbs += probs[i]
                
                newBaseProb = random.random()
                probSeen = 0
                for i in range(len(probs)):
                    if (i != consensusIndex):
                        probSeen += probs[i]/sumProbs
                        if(probSeen >= newBaseProb):
                            nextLetter = hmm.alphabet[i]
                            break
            modSeq += nextLetter
        return modSeq
            
    
    # places repeat in the sequence
    def insertRepeatSequence(self, repeat, sequence):
        return sequence + repeat
    
    # Create Repeat Entry to store in PSM list and same to file
    def createRepeatForPSM(self, start, end, HMM, modifiedSequence):
        psmRepeat = RptMatrixMod.RptMatrix(None,None)
        countMat = self.createCountMatrix(HMM, modifiedSequence )
        psmRepeat.M = countMat
        psmRepeat.start = start
        psmRepeat.finish = end
        psmRepeat.class_name = HMM.name.strip()
        return psmRepeat
        
    
    #Writes the PSM file to disk
    def writePSMFile(self, fileName, repeatList):
        pickle.dump(repeatList, open(fileName, "wb"))
    
    #Write the fasta sequence to file
    def writeFastaSequence(self,fileName, seq):
        outFile = open(fileName,"w")
        outFile.write(">" + fileName +" \n")
        outFile.write(seq)
    
    #Makes the count matrix for the repeat object
    def createCountMatrix(self,hmm, sequence):
        M = numpy.zeros((4,4))
        consensusSequence = self.getConsensusSequence(hmm)
        base2index = {'A':0,'C':1,'G':2,'T':3}                 # Dictionary mapping bases to their index position. 
        baseSet = {'A','C','G','T'} 
        for c1,c2 in zip(consensusSequence, sequence):
            if {c1,c2} <= baseSet:
                M[base2index[c1],base2index[c2]] += 1
        return M
    
    def getConsensusSequence(self, hmm):
        consensusSeq = ""
        for x in range(1,len(hmm.emissionProbs)):
            consensusSeq += (hmm.alphabet[hmm.getConsensusIndex(x)])
        return consensusSeq
    
    #Generate Random Sequence
    def generateRandomSequence(self, numChars = 1000):
        DNA=""
        for count in range(numChars):
            DNA+=random.choice("CGTA")
        return DNA
    
    #Create Fasta and PSM Files
    def populateFiles(self):
        numRegions = int((self.endPercentChange - self.startPercentChange)/self.percentIncrement + 1)
        percentChange = self.startPercentChange
        indexFile = open(self.repeatIndicesFile,"w")
        fastaSeq = ""
        repeatList = []
        for i in range(numRegions):
            for k in range(17):
                #Place Random Sequnce and repeat to find
                seq = self.generateRandomSequence()
                fastaSeq = self.insertRepeatSequence(seq, fastaSeq)
                begin = len(fastaSeq)
                seq = self.modifySequence(self.repeatToFind, percentChange)
                fastaSeq = self.insertRepeatSequence(seq, fastaSeq)
                end = len(fastaSeq)
                #print("RepeatNum: " + str(i*17+k) + ", begin: " + str(begin + 1) + ", end: " + str(end))
                indexFile.write(str(begin +1) + ", " + str(end) +"\n")
                #print(seq)
                #place in other repeat
                for j in self.repeatHMMs:
                    seq = self.generateRandomSequence()
                    fastaSeq = self.insertRepeatSequence(seq, fastaSeq)
                    repeatStart = len(fastaSeq)
                    seq = self.modifySequence(j, percentChange)
                    fastaSeq = self.insertRepeatSequence(seq, fastaSeq)
                    repeatEnd = len(fastaSeq)
                    repeatList.append(self.createRepeatForPSM(repeatStart,repeatEnd,j,seq))
            
            percentChange += self.percentIncrement
            
        self.writePSMFile("./PSMs/simulation.psm", repeatList)
        self.writeFastaSequence("./FastaFiles/simulation.fa", fastaSeq)
        
sc = SimulationContainer()          
sc.populateFiles()
                
                
    