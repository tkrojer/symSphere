import math
import glob
import os,sys

def getAlphabet():
    alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    return alphabet

def compareDict(refCAdict,symTMPdict,refList,symList,pdbname,cutoff):
    for a in refCAdict:
        x1 = refCAdict[a][0]
        y1 = refCAdict[a][1]
        z1 = refCAdict[a][2]
        for b in symTMPdict:
            x2 = symTMPdict[b][0]
            y2 = symTMPdict[b][1]
            z2 = symTMPdict[b][2]
            distance=math.sqrt(math.pow(float(x1)-float(x2),2)+math.pow(float(y1)-float(y2),2)+math.pow(float(z1)-float(z2),2))
            if distance < cutoff:
                if pdbname+b not in symList:
                    symList.append(pdbname+'#'+b)
                if a not in refList:
                    refList.append(a)
    return refList,symList



def getDict(pdb,which):
    dictTMP = {}
    for line in open(pdb):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atomname = str(line[13:16]).replace(' ','')
            resname=str(line[17:20]).replace(' ','')
            chainID=str(line[21:23]).replace(' ','')
            resseq=str(line[23:26]).replace(' ','')
            altLoc=str(line[16:17]).replace(' ','')
            occupancy=str(line[56:60]).replace(' ','')
            X=float(line[30:38])
            Y=float(line[38:46])
            Z=float(line[46:54])
            if which == 'CA':
                if atomname == 'CA':
                    if chainID+'#'+resseq not in dictTMP:
                        dictTMP[chainID+'#'+resseq] = [X,Y,Z]
            if which == 'all':
                if chainID+'#'+resseq+'#'+atomname not in dictTMP:
                    dictTMP[chainID+'#'+resseq+'#'+atomname] = [X,Y,Z]
    return dictTMP

def writeTMPpdb(pdbin,caList,type,newChain):
    out = ''
    pdbname = pdbin.replace('.pdb','')
    if type == 'ref' or type == 'sym':
        pre = 'tmp'
    else:
        pre = 'out'

    n = 0
    for line in open(pdbin):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            chainID=str(line[21:23]).replace(' ','')
            resseq=str(line[23:26]).replace(' ','')
            resname=str(line[17:20]).replace(' ','')
            atomname = str(line[13:16]).replace(' ','')
            if type == 'ref':
                if chainID+'#'+resseq in caList:
                    out += line
            if type == 'sym':
                if pdbname+'#'+chainID+'#'+resseq in caList:
                    out += line
            if type == 'other':
                if pdbname+'#'+chainID+'#'+resseq+'#'+atomname in caList:
                    l = list(line)
                    l[13] = 'O'
                    l[14] = ' '
                    l[15] = ' '
                    l[17] = 'H'
                    l[18] = 'O'
                    l[19] = 'H'
                    l[21] = newChain
                    l[77] = 'O'
                    if len(str(n+1)) == 1:
                        l[23] = ' '
                        l[24] = ' '
                        l[25] = str(n+1)
                    elif len(str(n+1)) == 2:
                        l[23] = ' '
                        l[24] = str(n+1)[0]
                        l[25] = str(n+1)[1]
                    elif len(str(n+1)) == 3:
                        l[23] = str(n+1)[0]
                        l[24] = str(n+1)[1]
                        l[25] = str(n+1)[2]
                    line = "".join(l)
                    n += 1
                    out += line


    f = open(pre+'_'+pdbin,'w')
    f.write(out)
    f.close()


def cleanup(symROOT):
    os.system('/bin/rm tmp_'+symROOT+'*.pdb')
    os.system('/bin/rm out_tmp_'+symROOT+'*.pdb')


if __name__ == '__main__':
    alpha = getAlphabet()

    ref = sys.argv[1]
    symROOT = sys.argv[2]
    cutoff = int(sys.argv[3])

    refCAdict = getDict(ref,'CA')
    refList = []
    symList = []
    nSym = len(glob.glob(symROOT+'*.pdb'))
    for n,symfile in enumerate(sorted(glob.glob(symROOT+'*.pdb'))):
        print 'reading ',n+1,' of ',nSym,' PDB files'
        pdbname = symfile.replace('.pdb','')
        symTMPdict = getDict(symfile,'CA')
        refList,symList = compareDict(refCAdict,symTMPdict,refList,symList,pdbname,cutoff+5)
        writeTMPpdb(symfile,symList,'sym',alpha[n])
    writeTMPpdb(ref,refList,'ref',alpha)
    refALLdict = getDict('tmp_'+ref,'all')
    for n,symfile in enumerate(sorted(glob.glob('tmp_'+symROOT+'*.pdb'))):
        print 'reading ',n+1,' of ',nSym,' PDB files'
        pdbname = symfile.replace('.pdb','')
        symTMPdict = getDict(symfile,'all')
        refList,symList = compareDict(refALLdict,symTMPdict,refList,symList,pdbname,cutoff)
        writeTMPpdb(symfile,symList,'other',alpha[n])
    out = ''
    for n,symfile in enumerate(sorted(glob.glob('out_tmp_'+symROOT+'*.pdb'))):
        for line in open(symfile):
            out +=line
    f = open('symAtoms.pdb','w')
    f.write(out)
    f.close()

    cleanup(symROOT)
