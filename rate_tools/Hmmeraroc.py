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
                    default = "./FastaFiles/chr22.fa", help='Set Fasta File')
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
    index = 0
    for part in fInfo.partitions:
        hmmFile = fInfo.hmmFolder + fInfo.hmmStub + str(index) 
        seqFile = fInfo.seqFolder + fInfo.fastaStub + str(index)
        resFile = results.resultsFolder + str(index) + ".res"
        hmmerInter = HmmerInterface.HmmerInfo(seqFile, hmmFile, resFile)
        hmmerInter.findResults()
        index +=1
        FileMaker.parseResultsFile(resFile, resFile + "cleaned", part.startIndex)
    