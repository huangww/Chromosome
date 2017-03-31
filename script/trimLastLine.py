import os
import glob


def trimFile(fname):
    """ trim the last line of a file

    :fname: the file name to treat
    :returns: None

    """
    f = open(fname, 'r')
    lineList = f.readlines()
    f.close()
    firstLine = lineList[0].split('\t')
    lastLine = lineList[-1].split('\t')
    if len(lastLine) == len(firstLine):
        pass
    else:
        f = open(fname, 'w')
        f.writelines([line for line in lineList[:-1]])
        f.close()
    return

def main():
    """ Main function of the script
    :returns: TODO

    """
    filePattern = 'BeadRod_rd_N100_F1_T*_*.dat'
    fileList = [fname for fname in glob.iglob(filePattern) if os.path.getsize(fname) > 1e5]
    for fname in fileList:
        trimFile(fname)
    return

if __name__ == "__main__":
    main()
