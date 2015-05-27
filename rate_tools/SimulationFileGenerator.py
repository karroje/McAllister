import HMMEditor
class SimulationContainer:
    def __init__(self, start = 0, end = .5,incr = .01):
        self.otherRepeatFileNames = ["./HMMs/LTR16a.hmm","./HMMs/MARNA.hmm","./HMMs/MER115.hmm","./HMMs/mlt1l.hmm", "./HMMs/tigger8.hmm"]
        self.repeatHMMs = self.retrieveRepeatList(self.otherRepeatFileNames)
        self.startPercentChange = start
        self.endPercentChange = end
        self.percentIncrement = incr
        self.repeatToFind = self.retrieveHMM("./HMMs/MIR.hmm")
        
        
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
        pass
    
    # places repeat in the sequence
    def insertRepeatSequence(self, repeat, sequence):
        return sequnece + repeat
    
    # Create Repeat Entry to store in PSM list and same to file
    def createRepeatForPSM(self, start, end, HMM, modifiedSequence):
        countMat - createCountMatrix(HMM, )
    
    #Writes the PSM file to disk
    def writePSMFile(self, fileName, repeatList):
        pickle.dump(repeatList, open(fileName, "wb"))
    
    #Write the fasta sequence to file
    def writeFastaSequence(self,fileName, seq):
        pass
    
    #Makes the count matrix for the repeat object
    def createCountMatrix(self,hmm, sequence):
        pass
    #Generate Random Sequence
    def generateRandomSequence(self, numChars = 1000):
        pass
    
    #Create Fasta and PSM Files
    def populateFiles(self):
        numRegions = (self.endPercentChange - self.startPercentChange)/self.percentIncrement + 1
        percentChange = start
        fastaSeq = ""
        repeatList = []
        for i in range(numRegions):
            for k in range(17):
                for j in self.repeatHMMs:
                    seq = self.generateRandomSequence()
                    fastaSeq = self.insertRepeatSequence(seq, fastaSeq)
                    repeatStart = len(fastaSeq)
                    seq = self.modifySequence(j, percentChange)
                    fastaSeq = self.insertRepeatSequence(seq, fastaSeq)
                    repeatEnd = len(fastaSeq)
                    repeatList.append(createRepeatForPSM(repeatStart,repeatEnd,j,seq))
            percentChange += self.percentIncrement
            
        self.writePSMFile("./PSMs/simulation.psm", repeatList)
        self.writeFastaSequence("./FastaFiles/simulation.fa", seq)
                
                
                
    