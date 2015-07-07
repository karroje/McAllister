import requests
import gzip
import StringIO
#import HMMEditor

def getRequestData(entry):
    data= {
        'file': 'hmm'
        }
    r = requests.post('http://www.dfam.org/download/model/' + entry, data = data)
    return r

def createTempFile(req, fileName):
    tempFile = open(fileName,'wb')
    #data = StringIO.StringIO(r.text)
    tempFile.write(req.content)
    tempFile.close()

def downloadFiles(startIndex,EndIndex):
    for i in range(startIndex,EndIndex):
        print i
        numWithZeros = "{0:04d}".format(i)
        req = getRequestData("DF000" + numWithZeros)
        createTempFile(req,'tempFile.hmm.gz')
        fh = unzipFile('tempFile.hmm.gz')
        saveNamedFile(fh)

def unzipFile(fileName):
    f = gzip.open('tempFile.hmm.gz', 'rb')
    return f

def saveNamedFile(fh):
    try:
        line1 = fh.readline()
        line2 = fh.readline()
        line2List = line2.split()
        name = line2List[1]
        hmmFile = open("./AllHMMs/" + name + ".hmm",'w')
        hmmFile.write(line1)
        hmmFile.write(line2)
        hmmFile.write(fh.read())
        #print file_content
        fh.close()
        hmmFile.close()
    except IOError:
        print("Error")

downloadFiles(1, 1338)
