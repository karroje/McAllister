import HmmerInterface
import FileMaker
import argparse
import datetime
import os.path

class FileInfo:
    def __init__(self):
        self.partitions = []
        self.folder = ""
        self.hmmFolder =""
        self.seqFolder =""
        self.hmmStub = ""
        self.fastaStub = ""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action='store', dest='fastaFile',
                    default = "./FastaFiles/chr11.fa", help='Set Fasta File')
    parser.add_argument('-m', action='store', dest='hmmFile',
                    default = "./HMMs/MIR.hmm", help='Set Orig HMM file')
    parser.add_argument('-p', action='store', dest='HMMERPATH',
                    default = "/usr/local/bin/", help='Path to HMM software')
    res = "./HMMResults" + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') +"/"
    parser.add_argument('-r', action='store', dest='resultsFolder',
                    default = res, help='Where results will be located')
    results = parser.parse_args()
    
    fInfo = FileInfo()
    FileMaker.createFiles(results.fastaFile, results.hmmFile, results.resultsFolder, fInfo)
    
    procs = []
    for index in xrange(len(fInfo.partitions)):
        part = fInfo.partitions[index]
        hmmFile = fInfo.hmmFolder + fInfo.hmmStub + str(index) 
        seqFile = fInfo.seqFolder + fInfo.fastaStub + str(index)
        resFile = results.resultsFolder + str(index) + ".res"
        hmmerInter = HmmerInterface.HmmerInfo(seqFile, hmmFile, resFile)
        procs.append(hmmerInter.findResults())
        #Run Original Hmm file
        hmmerInterOrig = HmmerInterface.HmmerInfo(seqFile, results.hmmFile, resFile + "orig")
        procs.append(hmmerInterOrig.findResults())
        procs[-2].wait()
        procs[-1].wait()
        
    repBasesMod = 0
    repBasesOrig = 0
    for index in xrange(len(fInfo.partitions)):
        p = procs[index*2]
        part = fInfo.partitions[index]
        resFile = results.resultsFolder + str(index) + ".res"
        #p.wait()
        modResults = FileMaker.parseResultsFile(resFile, resFile + "cleaned", part.startIndex)
        pOrig = procs[index*2 + 1]
        resFileOrig = results.resultsFolder + str(index) + ".resorig"
        #pOrig.wait()
        origResults =FileMaker.parseResultsFile(resFileOrig, resFileOrig + "cleaned", part.startIndex)
        repBasesMod += modResults[0]
        repBasesOrig += origResults[0]
        if(float(modResults[0])/modResults[1] > float(origResults[0])/origResults[1]):
            print( "More results in modified.  Partition: " + str(index))
        elif(float(modResults[0])/modResults[1] < float(origResults[0])/origResults[1]):
            print("Less results in modified.  Partition: " + str(index))
        else:
            print("Same results.  Partition: " + str(index))
        print(str(modResults[0]) + " : " + str(origResults[0]) + " : "  + str(part.calculateChange()))
    print(str(repBasesMod) + " : " + str(repBasesOrig))
        
    