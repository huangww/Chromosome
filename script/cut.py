import os
import glob

beadNumber = 98

dataDir = 'data/'
# fileString = 'r_N98_*.dat'
# fileList = glob.glob(dataDir + fileString)

fileName = fileList[0]

newfile = open(fileName[:-4]+'_cut.dat', 'w')

lineCount = 0
commentLinePre = 0
frame = []
with open(fileName) as f:
    for line in f:
        if line[0] == '#':
            commentLineNow = lineCount
            frameLen = commentLineNow - commentLinePre
            if frameLen == beadNumber+1:
                newfile.writelines(frame)
            commentLinePre = lineCount
            frame[:] = []
        frame.append(line)
        lineCount += 1

frameLen = lineCount - commentLinePre
if frameLen == beadNumber+1:
    newfile.writelines(frame)

newfile.close()

