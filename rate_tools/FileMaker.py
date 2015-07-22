"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import uuid
import os
import pickle
import copy
import DNAPartitioner
import HMMEditor
import RptMatrixMod
import RMfileReader




charRead = 0

def maketempDirectory():
    tempFileName = os.path.join(os.path.dirname(__file__))
    tempFileName += '/' + str(uuid.uuid4()) +'/'
    #print(tempFileName)
    os.makedirs(tempFileName)
    return tempFileName

def aggregrateResultFiles(resultsFolder, ext):
    aggFile = open(resultsFolder+"Cumulative" +ext,'w')
    fileNames = os.listdir(resultsFolder)
    fileNum = 0
    fileCheck = str(fileNum) + ext
    while ( fileCheck in fileNames):
        fh = open(resultsFolder+fileCheck,'r')
        lineNum = 0
        for line in fh:
            if lineNum > 2 and line.strip():
                aggFile.write(line)
            lineNum +=1
        
        fileNum +=1
        fileCheck = str(fileNum) + ext
         
    
def parseResultsFile(infile, outFile, startIndex):
    summary = False
    startLine = "Scores for complete hits:"
    endLine1 = "------ inclusion threshold ------"
    endLine2 = "Annotation for"
    fh = open(infile, 'r')
    fw = open(outFile,'w')
    totalBases = 0.0
    repeatBases = 0.0
    for line in fh:
        if(startLine in line):
            summary=True
        elif(endLine1 in line or endLine2 in line):
            summary=False
        
        if (summary):
            repeatBases += countBases(line)
            line = changeBases(line, startIndex)
            fw.write(line)
        if("Target sequences" in line):
            totalBases = extractTotalBases(line)
    return [repeatBases,totalBases]
            
def extractTotalBases(line):
    totalBases = -1
    pieces = line.split()
    if(len(pieces)>3):
        end = pieces[3][1:]
        endPieces = end.split()
        if(len(endPieces)>0):
            if(isInt(endPieces[0])):
                totalBases = int(endPieces[0])
    return totalBases
    
            
def changeBases(line, offset):
    pieces = line.split()
    if(validResultLine(line)):
        pieces[5] = str(int(pieces[5]) + offset)
        pieces[4] = str(int(pieces[4]) + offset)
    line = "   ".join(pieces) + "\n"
    return line

def countBases(line):
    bases = 0
    if(validResultLine(line)):
        result = line.split()
        bases = abs(int(result[5]) - int(result[4])) + 1
    return bases

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
def validResultLine(line):
    valid = True
    result = line.split()
    if(len(result) < 6):
        valid = False
    elif(not isInt(result[4])):
        valid = False
    elif(not isInt(result[5])):
        valid = False
    return valid

def makeFastaFile(part, fh, fileName, lastCharIndex,seqFolder):
    distFromPrev = max(part.startIndex - lastCharIndex,0)
    fh.read(distFromPrev)
    global charRead 
    charRead += distFromPrev
    #print(distFromPrev)
    fh2 = open(seqFolder+fileName,'w')
    fh2.write("> " + fileName + "from " + str(part.startIndex)
              + " to " + str(part.endIndex) + "\n")
    s = fh.read(part.endIndex - part.startIndex + 1)
    #print(len(s))
    fh2.write(s)
    fh2.close
    charRead += part.endIndex - part.startIndex + 1
    #print(charRead)
    return part.endIndex + 1
    

def populateFastaFiles(partitionList,seqFolder, startFile):
    fh = open(startFile,'r')
    fh.readline()
    index = 0
    lastCharIndex = 0
    lastPart = startFile.split("/")[-1].split(".")[0]
    for part in partitionList:
        lastCharIndex = makeFastaFile(part, fh, lastPart + "Part" + str(index),
                                       lastCharIndex, seqFolder)
        index +=1
        #print(index)
    fh.close()
    return lastPart + "Part" 
       
def populateHMMFiles(partitions, hmmFolder, origHmm, startFile, negOnly =False):
    index = 0
    lastPart = startFile.split("/")[-1].split(".")[0]
    for part in partitions:
        newHmm = copy.deepcopy(origHmm)
        #[:-1] is used to remove newline character from name
        newHmm.name = newHmm.name[:-1] + lastPart + "Part" + str(index)
        changeProb = part.calculateChange()
        if(negOnly):
            changeProb = min(0, changeProb)
        newHmm.addToDecimalEmissionProbs(changeProb)
        newHmm.convertDecimalProbsToEmission()
        newHmm.writeFile(hmmFolder + newHmm.name)
        index +=1
    return  origHmm.name[:-1] + lastPart + "Part" 
    
def createFiles(startingFasta, startingHMM, psmFile, resFolder=None, fileInfo = None,
                negOnly = False):
    folder = maketempDirectory()
    hmmFolder = folder+"Hmms/"
    os.makedirs(hmmFolder)
    seqFolder = folder+"Sequences/"
    os.makedirs(seqFolder)
    if (resFolder != None):
        os.makedirs(resFolder)
    hmmName = startingHMM.split("/")[-1].split(".")[0]
    partitioner = DNAPartitioner.PartitionMaker(search = hmmName, psm=psmFile)
    partitioner.createPartitions() 
    partitioner.printPartitionList()
    partFile = open("partFile", "wb")
    #pickle.dump(partitioner,partFile)
    #partFileRead =  open("partFile",'r')
    #partitioner = pickle.load(partFileRead)
    hmmOrig = HMMEditor.hmmprofile()
    hmmOrig.parseFile(startingHMM)
    fastaStub = populateFastaFiles(partitioner.partitionList, seqFolder, startingFasta)
    hmmStub = populateHMMFiles(partitioner.partitionList, hmmFolder, hmmOrig, startingFasta, negOnly)
    if(fileInfo != None):
        fileInfo.folder=folder
        fileInfo.partitions = partitioner.partitionList
        fileInfo.hmmFolder = hmmFolder
        fileInfo.seqFolder = seqFolder
        fileInfo.fastaStub = fastaStub
        fileInfo.hmmStub = hmmStub
    
    
def createPickledFiles():
    RptMatrixMod.create_psm("/Users/MikeMcAllister/cache/human/hg18/seq/rmsk")
#print(createPickledFiles())
if __name__ == "__main__":
    createPickledFiles()
    