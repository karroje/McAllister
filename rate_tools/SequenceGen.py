"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import random
def weightedchoice(items): # this doesn't require the numbers to add up to 100
    return random.choice("".join(x * y for x, y in items))

def changeChars(perChange, strToChange):
    newString =""
    for c in strToChange:
        if random.random() < perChange:
            choices = [("A",1),("C",1),("G",1),("T",1)]
            choices.remove((c,1))
            c = weightedchoice(choices)
        newString +=c
    return newString
            
def makeFastaSequence(outputFile = "testDNASeq", repeatSeq ="TTGCAATACACAAGTGATCG", 
                      len = 500, numRepeats=25, maxChange =.5): 
    outFile = open(outputFile,"w")
    random.seed(1)
    repeats = 0
    outFile.write(">" + outputFile +" \n")
    for i in range(0,len):
        if i % (len/numRepeats) == 0:
            outString = changeChars(maxChange*repeats/numRepeats,repeatSeq)
            outFile.write(outString)
            repeats +=1
            outFile.write("\n");
        else:
            choices = [("A",1),("C",1),("G",1),("T",1)]
            c = weightedchoice(choices)      
            outFile.write(c)
def makeAlignmentSequences(outputFile = "testAlign", repeatSeq ="TTGCAATACACAAGTGATCG", 
                    numAligns=5, maxChange =.5): 
    outFile = open(outputFile,"w")
    random.seed(1)
    alignments = 0
    outFile.write("# STOCKHOLM 1.0 \n")
    for i in range(0,numAligns):
        outString = changeChars(maxChange*alignments/numAligns,repeatSeq)
        outFile.write("Sample" + str(alignments) + "     ")
        outFile.write(outString)
        alignments +=1
        outFile.write("\n");
    outFile.write("//")
makeAlignmentSequences()