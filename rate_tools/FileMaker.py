import uuid
import os
import DNAPartitioner
import RptMatrixMod
import RMfileReader
from RMfileReader import Repeat





def maketempDirectory():
    tempFileName = os.path.join(os.path.dirname(__file__))
    tempFileName += '/' + str(uuid.uuid4()) +'/'
    print(tempFileName)
    os.makedirs(tempFileName)
    return tempFileName



def parseHMMFile(infile, outFile, startIndex):
    summary = False
    startLine = "Scores for complete hits:"
    endLine1 = "------ inclusion threshold ------"
    endLine2 = "Annotation for"
    fh = open(infile, 'r')
    fw = open(outFile,'w')
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
    line = "\t".join(pieces) + "\n"
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
    
    
def createFiles():
    folder = maketempDirectory()
    os.makedirs(folder+"Hmms/")
    os.makedirs(folder+"Sequences/")
    partitioner = DNAPartitioner.PartitionMaker()
    partitioner.createPartitions() 
    return(parseHMMFile("results", folder + "smallResults.txt", 2847))
    
def createPickledFiles():
    RptMatrixMod.create_psm("/Users/MikeMcAllister/cache/human/hg18/seq/rmsk")
#print(createPickledFiles())
print(createFiles())