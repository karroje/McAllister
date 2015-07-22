"""
Copyright (c) 2015, Michael McAllister and Dr. John Karro
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

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
