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