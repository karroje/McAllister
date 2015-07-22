"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import math
from numpy import zeros
from asymm_tools import compute_d
from asymm_tools import compute_P
    
class hmmprofile:
    def __init__(self):
        self.name = "tempName"
        self.repeatLength = 0
        self.firstline =""
        self.transitionLine = ""
        self.header = ""
        self.alphabet=[]
        #List of list the first 4 numbers are the transition probs the 6th is the consensus letter
        self.emissionProbs=[]
        #List of list one for each letter of alphabet in order
        self.insertionProbs=[]
        # List of transition probabilities see transitionLine for more details
        self.transitionProbs=[]
        self.emissionDecimalProbs=[]
    
    def defineAlphabet(self,line):
        splitLine = line.split()
        self.alphabet = splitLine[1:]
    
    def convertEmissionProbsToDecimal(self):
        #print self.emissionProbs
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
                self.emissionProbs[x][y] = str(round(negLN,5))
        #print self.emissionProbs
        
    def getConsensusIndex(self, indexNum):
        letter = self.emissionProbs[indexNum][5]
        letter = letter.upper()
        index = -1
        try:
            index = self.alphabet.index(letter)
        except Exception:
            pass
        return index
    
    def assignProportions(self, probs, consensusIndex, amtToAdd):
        sum = 0
        for i in range(len(probs)):
            if (i != consensusIndex):
                sum += probs[i]
        deltas = []
        for i in range(len(probs)):
            if (i != consensusIndex):
                deltas.append(-1*amtToAdd*probs[i]/sum)
            else:
                deltas.append(amtToAdd)
        #print "Original Deltas"
        #print deltas
        return deltas
    """
    Reallocates extra from probabilities below 0 or over 100.  Should not be used often
    hence the print statement
    """
    def reallocateExtras(self, uneditSet, amtToAdd, deltas):
        print"Reallocating percentages due to crossing bottom or top Threshold percentage"
        partitions = len(deltas) - len(uneditSet)
        for i in range(len(deltas)):
            if i not in uneditSet:
                deltas[i] += amtToAdd/partitions
        return deltas
    
    def calibrateDeltasToThresholds(self, conIndex, topThreshold, bottomThreshold, origProbs, startingDeltas
                                    , altsAtThreshold = None):
        changeToAlts = 0
        if (altsAtThreshold == None):
            altsAtThreshold = set()
        for x in range(len(origProbs)):
            origProb = origProbs[x]
            deltaProb = startingDeltas[x]
            newDeltaProb = deltaProb
            if origProb + deltaProb > topThreshold:
                newDeltaProb = topThreshold - origProb 
                altsAtThreshold.add(x)              
            elif origProb + deltaProb < bottomThreshold:
                newDeltaProb = bottomThreshold - origProb
                altsAtThreshold.add(x)
            changeToAlts += deltaProb - newDeltaProb
            startingDeltas[x] = newDeltaProb
        if changeToAlts > 0 or changeToAlts < 0:
            altsAtThreshold.add(conIndex)
            self.reallocateExtras(altsAtThreshold, changeToAlts, startingDeltas)
            self.calibrateDeltasToThresholds(conIndex, topThreshold, bottomThreshold, origProbs, startingDeltas, altsAtThreshold)
        return startingDeltas
    '''
    Ensures that the consensus match still has the highest probability by taking average of all elements 
    higher than it and setting each element equal to the average
    '''
    def maintainConsensus(self, deltas, origProbs, consensusIndex):
        newConsensusProb = origProbs[consensusIndex] + deltas[consensusIndex]
        sumOfHigherProbs = 0
        indexesOfHigherProbs = []
        for i in range(len(origProbs)):
            newProb = origProbs[i] + deltas[i]
            if newProb > newConsensusProb:
                sumOfHigherProbs += newProb
                indexesOfHigherProbs.append(i)
        if sumOfHigherProbs > 0:
            sumOfHigherProbs +=newConsensusProb
            indexesOfHigherProbs.append(consensusIndex)
            avgProb = sumOfHigherProbs/len(indexesOfHigherProbs)
            for i in indexesOfHigherProbs:
                deltas[i] = avgProb - origProbs[i]
        return deltas
    
    '''
    Determines if any of the following are violated and shifts probabilites accordingly:
    1)something goes over 100% ... This is probably an error, but shifts to other alternatives
    2) Something goes under 0% ... This is probably okay, but shifts to other alternatives
    3) A letter other than the consensus has the highest percentile ... assignes each the average between
    the set of letters where this would be applicable
    '''
    def determineDeltaChanges(self, origProbs, consensusIndex, amtToAdd, proportional, topThreshold, bottomThreshold):
        startingDeltas = []
        if proportional:
            startingDeltas = self.assignProportions(origProbs, consensusIndex, amtToAdd)
        else:
            startingDeltas = self.assignProportions([1]*len(origProbs), consensusIndex, amtToAdd)
        # make sure thresholds aren't violated and calibrate accordingly by moving alternatives
        deltas = self.calibrateDeltasToThresholds(consensusIndex, topThreshold, bottomThreshold, origProbs, startingDeltas)
        self.maintainConsensus(deltas, origProbs, consensusIndex)
        return deltas
            
           
            
            
    
    '''    
    the parameter amtToAdd is the amount that will be added to the most likely value
    To reduce the value just send in a negative number
    The parameter proportional indicates whether to reduce/increase the other proportions
    in relation to their existing value for or uniformly across all alternatives
    E.g. if A = 60%, C=20%, G =10%, T = 10% ... if you add -10% to A and have proportional you will
    add 5% to C and 2.5 to each of G and T.  If not proportional, then you will will add 3.33 to each
    
    Other notes the method checks that the initial value does not go to 100% nor does it drop lower
    than any other letter.   If it does go below any other letters each letter is set to the average 
    of the set.  These checks should be irrelevant in most common cases 
    '''
    def addToDecimalEmissionProbs(self, amtToAdd, topThreshold = .96, bottomThreshold = .01, proportional = True):
        #check to see if input is reasonable
        #print self.emissionDecimalProbs
        #print("amtToAdd: " + str(amtToAdd))
        if amtToAdd > .75 or amtToAdd < -.75:
            raise Exception("The amount of change is outside the expected range")
        if topThreshold > 1 or topThreshold < .26:
            raise Exception("The topThreshold is outside the expected range")
        if bottomThreshold > .25 or bottomThreshold < .00001:
            raise Exception("The bottomThreshold is outside the expected range")
        for i in range(1,len(self.emissionDecimalProbs)):
            index = self.getConsensusIndex(i)
            #print(i)
            if index > -1:
                currentProb = self.emissionDecimalProbs[i][index]
                amtToAddIndex = amtToAdd
                if currentProb + amtToAdd > topThreshold:
                    amtToAddIndex = topThreshold - currentProb
                #print(i)
                deltas = self.determineDeltaChanges(self.emissionDecimalProbs[i], index, amtToAddIndex,
                                                     proportional, topThreshold, bottomThreshold)
                for j in range(len(self.emissionDecimalProbs[i])):
                    self.emissionDecimalProbs[i][j] += deltas[j]
        #print self.emissionDecimalProbs
        
                
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
        self.repeatLength = len(self.emissionDecimalProbs) - 1
          
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
        file.close()
    def writeConsensusSequence(self, fileName):
        file = open(fileName,"w")
        file.write("> " + fileName + "\n")
        for x in range(1,len(self.emissionProbs)):
            base ='N'
            consensusIndex = self.getConsensusIndex(x)
            if consensusIndex > -1:
                base = self.alphabet[consensusIndex]
            file.write(base)
    
    def generateExpectedCountMatrix(self):
        cntMat = zeros((4,4))
        for x in range(1,len(self.emissionDecimalProbs)):
            row = self.getConsensusIndex(x)
            if row > -1:
                for j in range(4):
                    cntMat[row][j] = self.emissionDecimalProbs[x][j]*100
        return cntMat

def runTest():
    #testProfile = hmmprofile()
    #testProfile.parseFile("./HMMs/ERV24_Prim.hmm")
    #print(testProfile.emissionDecimalProbs)
    #print(testProfile.generateExpectedCountMatrix())
    #testProfile.addToDecimalEmissionProbs(-.7)
    #testProfile.addToDecimalEmissionProbs(.7)
    #testProfile.convertDecimalProbsToEmission()
    #testProfile.writeFile("HMMGreaterChange")
    testProfile2 = hmmprofile()
    testProfile2.parseFile("./AllHmms/MER90.hmm")
    testProfile2.addToDecimalEmissionProbs(0.0106380591917)
    #testProfile2.writeConsensusSequence("./HMMs/MADE2.fa")
    #print(testProfile2.emissionDecimalProbs)
    #print(testProfile2.generateExpectedCountMatrix())
if __name__ == "__main__":
    print(runTest())