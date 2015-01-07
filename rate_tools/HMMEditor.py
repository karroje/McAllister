import math

class hmmprofie:
    def __init__(self):
        self.name = "tempName"
        self.firstline =""
        self.transitionLine = ""
        self.header = ""
        self.alphabet=[]
        self.emissionProbs=[]
        self.insertionProbs=[]
        self.transitionProbs=[]
        self.emissionDecimalProbs=[]
    def defineAlphabet(self,line):
        splitLine = line.split()
        self.alphabet = splitLine[1:]
    
    def convertEmissionProbsToDecimal(self):
        print self.emissionProbs
        self.emissionDecimalProbs= []
        for x in range(len(self.emissionProbs)):
            self.emissionDecimalProbs.append([])
            for y in range(len(self.alphabet)):
                negLN = float(self.emissionProbs[x][y])
                decimal = math.exp(-1*negLN)
                self.emissionDecimalProbs[x].append(decimal)
        self.testDecimalProbs()
    
    def convertDecimalProbsToEmission(self):
        for x in range(len(self.emissionProbs)):
            for y in range(len(self.alphabet)):
                LN = math.log(self.emissionDecimalProbs[x][y])
                negLN = -1*LN
                self.emissionProbs[x][y] = str(negLN)
        print self.emissionProbs
                
    def readProbs(self, file):
        tmpLine = file.readline()
        while tmpLine and tmpLine !="//\n":
            self.emissionProbs.append(tmpLine.split()[1:])
            tmpLine = file.readline()
            self.insertionProbs.append(tmpLine.split()[:])
            tmpLine = file.readline()
            self.transitionProbs.append(tmpLine.split()[:])
            tmpLine = file.readline()
        self.convertEmissionProbsToDecimal()
        self.convertDecimalProbsToEmission()    
    def testDecimalProbs(self):
        for x in self.emissionDecimalProbs:
            sum = 0;
            for y in x:
                sum += y
            if sum <.99 or sum>1.01:
                raise Exception("not valid probabilites")
    
    def parseFile(self,fileName):
        file = open(fileName,"r")
        self.firstline = file.readline()
        self.name = file.readline().split(" ")[-1]
        tmpLine = file.readline()
        threeChars = tmpLine[0:3]
        while threeChars != "HMM":
            self.header += tmpLine
            tmpLine = file.readline()
            threeChars = tmpLine[0:3]
        self.defineAlphabet(tmpLine)
        self.transitionLine = file.readline()
        self.readProbs(file)
    def writeAlphabet(self,file):
        #ten spaces as per sample
        tmpStr="HMM          "
        for i in self.alphabet:
            tmpStr += i + "        "
        tmpStr += "\n"
        file.write(tmpStr)
        
    def writeProbs(self,file):
        # ten spaces plus first line header
        tmpStr = "  COMPO   " + "  ".join(self.emissionProbs[0]) + "\n"
        tmpStr += "          " + "  ".join(self.insertionProbs[0]) + "\n"
        tmpStr += "          " + "  ".join(self.transitionProbs[0]) + "\n"
        file.write(tmpStr)
        for i in range (1, len(self.emissionProbs)):
            tmpStr=""
            tmpStr += " " * (7-len(str(i))) + str(i) + "   "
            tmpStr += "  ".join(self.emissionProbs[i]) + "\n"
            tmpStr += "          " + "  ".join(self.insertionProbs[i]) + "\n"
            tmpStr += "          " + "  ".join(self.transitionProbs[i]) + "\n"
            file.write(tmpStr)
    
    def writeFile(self,fileName):
        file = open(fileName,"w")
        file.write(self.firstline)
        file.write("NAME  " + self.name)
        file.write(self.header)
        self.writeAlphabet(file)
        file.write(self.transitionLine)
        self.writeProbs(file)
        file.write("//\n")
    
testProfile = hmmprofie()
testProfile.parseFile("testHmm")
testProfile.writeFile("testHMMout")