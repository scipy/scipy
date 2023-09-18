# Copyright (c) 2012 Ecma International.  All rights reserved.
# This code is governed by the BSD license found in the LICENSE file.

#--Imports---------------------------------------------------------------------
import argparse
import os
import sys
import re
import codecs

#--Globals---------------------------------------------------------------------
ECMA_LICENSE = '''/// Copyright (c) 2012 Ecma International.  All rights reserved.
/// This code is governed by the BSD license found in the LICENSE file.
'''

NEW_LICENSE_FIRST_LINE = re.compile(r"Copyright\s+\(c\)\s+20[0-9][0-9]\s+Ecma\s+International")
OLD_LICENSE_FIRST_LINE = re.compile(r"(Copyright\s+20[0-9][0-9]\s+Google)|(the\s+Sputnik\s+authors)|(Microsoft\s+Corporation)")
OLD_LICENSE_LAST_LINE  = re.compile(r"(ADVISED\s+OF\s+THE\s+POSSIBILITY\s+OF\s+SUCH\s+DAMAGE)|(This\s+code\s+is\s+governed\s+by\s+the\s+BSD\s+license\s+found\s+in\s+the\s+LICENSE\s+file)")

#Dirty way of determining if the contribution stems from Google or Microsoft
GOOGLE_RE = re.compile(r"[\\/]S([0-9]+)|(bp)(\.|_)[^\\/]+\.js$")
IETC_RE    = re.compile(r"[\\/][0-9]+(\.|_)[^\\/]+\.js$")

DEBUG = False
#------------------------------------------------------------------------------
def getAllJSFiles(dirName):
    '''
    Returns all JS files under dirName
    '''
    retVal = []
    if os.path.isfile(dirName) and dirName.endswith(".js"):
        retVal = [dirName]
    elif os.path.isdir(dirName):
        tempList = [os.path.join(dirName, x) for x in os.listdir(dirName)]
        for x in tempList:
            retVal += getAllJSFiles(x)
    #else:
    #    raise Exception("getAllJSFiles: encountered a non-file/non-dir:" + dirName)
    return retVal

#------------------------------------------------------------------------------
def handleFile(filePath):
    '''
    '''
    with open(filePath, "rb") as f:
        origLines = f.readlines()

    #See if it's already there
    if NEW_LICENSE_FIRST_LINE.search(origLines[0])!=None:
        #print "\talready there:\t", filePath
        return
    #TODO: Google employee needs to remove this elif
    #      and fix the next elif clause
    elif GOOGLE_RE.search(filePath)!=None:
        if DEBUG:
            print "\tignoring Google sources:\t", filePath
        return
    elif (IETC_RE.search(filePath))==None: #and (GOOGLE_RE.search(filePath)==None):
        errMsg = "\tno idea which license should be used for:\t" + filePath
        raise Exception(errMsg)
        return

    with codecs.open(filePath,'r','utf8') as f:
        bomPresent = f.read(2).startswith(u"\ufeff")
        if bomPresent:
            print "\tnon-ASCII file detected. Please modify by hand:", filePath
            return

    with open(filePath, "wb") as f:
        if DEBUG:
            print "\tmodified:\t", filePath
        #TODO: this isn't good enough...
        #if bomPresent:
        #    print "\tBOM was detected for:", filePath
        #    f.write(u"\ufeff")
        f.write(ECMA_LICENSE)

        writeIt = False
        for line in origLines:
            if writeIt:
                f.write(line)
            elif OLD_LICENSE_LAST_LINE.search(line)!=None:
                writeIt = True

        if not writeIt:
            print "\tError - didn't find end of the original license:\t", filePath

#--Main------------------------------------------------------------------------
if __name__=="__main__":
    __parser = argparse.ArgumentParser(description='Tool used to fix test file license headers')
    __parser.add_argument('tpath', action='store',
                          help='Full path to test cases. E.g., C:\repos\test262-msft\test\suite')
    ARGS = __parser.parse_args()
    if not os.path.exists(ARGS.tpath):
        print "Cannot fix tests in '%s' when it doesn't exist!" % ARGS.tpath
        sys.exit(1)

    ALL_JS_FILES = getAllJSFiles(ARGS.tpath)
    for fileName in ALL_JS_FILES:
        handleFile(fileName)
    print "Done!"
