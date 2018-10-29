#!/usr/bin/python
import copy


NjNbDict = {'2': ['0','1','2'],
            '3-4': ['0','1','2','>=3'],
            '5-6': ['0','1','2','>=3'],
            '7-8': ['0','1','2','>=3'],
            '9+': ['0','1','2','>=3'],
}

Nj = ['2','3-4','5-6','7-8','9+']
Nb = ['0','1','2','>=3']

f = open('../datFiles/DY_signal.dat', 'r')

rowList = []

for line in f:
    rowList.append(line.split())

## get rid of header
rowList.pop(0)

## initialize Muon yield dictionary 
NmuYield = {}

## setup dictionary format
for row in rowList:
    NmuYield[row[0]+":"+row[2]] = []

## deep copy muon dictonary format as other table elements
NelYield = copy.deepcopy(NmuYield); Extrap = copy.deepcopy(NmuYield)
Stat = copy.deepcopy(NmuYield); MCstat = copy.deepcopy(NmuYield) 
SysUp = copy.deepcopy(NmuYield); SysDn = copy.deepcopy(NmuYield); 
SysKin = copy.deepcopy(NmuYield); SysPur = copy.deepcopy(NmuYield)
TTZ = copy.deepcopy(NmuYield)
for row in rowList:
    NmuYield[row[0]+":"+row[2]].append(int(float(row[4])))
    NelYield[row[0]+":"+row[2]].append(int(float(row[6])))
    Extrap[row[0]+":"+row[2]].append(float(row[8]))
    Stat[row[0]+":"+row[2]].append(float(row[10]))
    MCstat[row[0]+":"+row[2]].append(float(row[12]))
    TTZ[row[0]+":"+row[2]].append(float(row[14]))
    SysUp[row[0]+":"+row[2]].append(float(row[16]))
    SysDn[row[0]+":"+row[2]].append(float(row[18]))
    SysKin[row[0]+":"+row[2]].append(float(row[20]))
    SysPur[row[0]+":"+row[2]].append(float(row[22]))

Bin = 0

for nj in Nj:
    print "\hline"
    for nb in Nb:
        if(nj+nb == '2>=3'):
            continue
        Bin+=1
        NjNbBin = nj+':'+nb
        # get total or average values for each Nj X Nb bin
        nmu = str(sum(NmuYield[NjNbBin]))
        nel = str(sum(NelYield[NjNbBin]))
        extrap = str(round(sum(Extrap[NjNbBin])/len(Extrap[NjNbBin]),3))
        stat = str(int(round(sum(Stat[NjNbBin])/len(Stat[NjNbBin]),2)*100))
        mcstat = str(int(round(sum(MCstat[NjNbBin])/len(MCstat[NjNbBin]),3)*100))
        sysup = str(int(round(sum(SysUp[NjNbBin])/len(SysUp[NjNbBin]),3)*100))
        sysdn = str(int(round(sum(SysDn[NjNbBin])/len(SysDn[NjNbBin]),3)*100))
        syskin = str(int(sum(SysKin[NjNbBin])/len(SysKin[NjNbBin])*100))
        syspur = str(int(round(sum(SysPur[NjNbBin])/len(SysPur[NjNbBin]),3)*100))
        ttz = str(int(round(sum(TTZ[NjNbBin])/len(TTZ[NjNbBin]),3)*100))
        print str(Bin)+" & "+nmu+" & "+nel+" & "+extrap+" & "+stat+" & "+syspur+" & $\pm  "+mcstat+"^{ +"+sysup+"}_{ -"+sysdn+"}$ & "+ttz+" & "+syskin+ " \\\\ "


f.close()
