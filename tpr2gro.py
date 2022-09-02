import numpy as np
import pandas as pd
import re
import time
import sys

startT=time.time()
# 1.read (gmx dump -s md.tpr > md.out ) file
with open(sys.argv[1]) as f:
    lines=f.readlines()

print("********************************")
print("input file is {}".format(sys.argv[1]))


# 2. find all molblocks
molblockstart=[]
box=[]
for line in lines:
    if line.startswith("   molblock"):
        molblockstart.append(lines.index(line))
    if line.startswith("box (3x3):"):
        box.append(re.findall(re.compile(r"=\{([\,\s\+\-e\.0-9]+)\}"),lines[lines.index(line)+1])[0].split(","))
        box.append(re.findall(re.compile(r"=\{([\,\s\+\-e\.0-9]+)\}"),lines[lines.index(line)+2])[0].split(","))
        box.append(re.findall(re.compile(r"=\{([\,\s\+\-e\.0-9]+)\}"),lines[lines.index(line)+3])[0].split(","))
totalatom={}
for i,j in enumerate(molblockstart):
    totalatom[re.findall(re.compile(r"\s[0-9]+\s\"([0-9a-zA-Z\-\_\+]+)\""),lines[j+1])[0]]=re.findall(re.compile(r"\s=\s([0-9]+)"),lines[j+2])[0]
    
    
# 3.read coord
for line in lines:
    if line.startswith("x ("):
        xlen=re.findall(re.compile(r"\s\(([0-9]+)x"),line)
        xstart=lines.index(line)+1
        xend=int(xlen[0])+xstart
coord=[re.findall(re.compile(r"\{([\s0-9a-zA-Ze\.\+\-\,]+)\}"),line)[0].split(",") for line in lines[xstart:xend]]
coord=pd.DataFrame(coord,columns=["x","y","z"])

# 4. read type mass charge atomidx resid atomname atomtype 
# 修改resindex思路，根据循环
# 先设置一个初始residx_tmp=0,在第一层循环时即对蛋白序列都只加上0
# 第二次循环时，residx_tmp取存储在tmp中的最后一个残基index进行累加
# 对于sol及NA离子，通过在循环内添加countindx实现自增
tmp=[]
type=re.compile(r"type=([\s0-9]+)")
mass=re.compile(r"m=([\s\.0-9e\+\-]+)")
charge=re.compile(r"q=([\s\.0-9e\+\-]+)")
adx=re.compile(r"atom\[([\s0-9]+)")
resid=re.compile(r"resind=([\s0-9\+\-]+)")
atn=re.compile(r"atomnumber=([\s0-9\+\-]+)")
atomname=re.compile(r"{name=\"([A-Za-z0-9\*\_]+)\"")
atomtype=re.compile(r"{name=\"([A-Za-z0-9\*\_]+)\",")
residuestart=[]
residx_tmp=0
for line in lines:
    for i,j in enumerate(totalatom):
        if line.startswith("      name="+"\""+j+"\""):
            Len=re.findall(re.compile(r"\s\(([0-9]+)"),lines[lines.index(line)+2])
            start=lines.index(line)+3
            end=start+int(Len[0])
            resn={}
            resns=re.findall(re.compile(r"\s\(([0-9]+)\)"),lines[start+3*int(Len[0])+2])[0]
            for res in range(int(resns)):
                resn[res]=re.findall(re.compile(r"{name=\"([a-zA-Z0-9\_\+\*]+)"),lines[start+3*int(Len[0])+res+3])[0]
            if int(totalatom[j])!=1:
                coutindx=1
                for n in range(int(totalatom[j])):
                    count=1
                    for m in lines[start:end]:
                        tmp.append([
                            re.findall(type,m)[0].split()[0],
                            re.findall(mass,m)[0].split()[0],
                            re.findall(charge,m)[0].split()[0],
                            re.findall(adx,m)[0].split()[0],
                            int(re.findall(resid,m)[0].split()[0])+coutindx+residx_tmp,
                            re.findall(atn,m)[0].split()[0],
                            re.findall(atomname,lines[start+int(Len[0])+count])[0].split()[0],
                            re.findall(atomtype,lines[start+2*int(Len[0])+count+1])[0].split()[0],
                            resn[int(re.findall(resid,m)[0].split()[0])],
                            j,
                        ])
                        count=count+1
                    coutindx=coutindx+1
                residx_tmp=tmp[-1][4]
            else:
                count=1
                for m in lines[start:end]:
                        tmp.append([
                            re.findall(type,m)[0].split()[0],
                            re.findall(mass,m)[0].split()[0],
                            re.findall(charge,m)[0].split()[0],
                            re.findall(adx,m)[0].split()[0],
                            int(re.findall(resid,m)[0].split()[0])+1+residx_tmp,
                            re.findall(atn,m)[0].split()[0],
                            re.findall(atomname,lines[start+int(Len[0])+count])[0].split()[0],
                            re.findall(atomtype,lines[start+2*int(Len[0])+count+1])[0].split()[0],
                            resn[int(re.findall(resid,m)[0].split()[0])],
                            j
                        ])
                        count=count+1
                residx_tmp=tmp[-1][4]
info=pd.DataFrame(tmp,columns=["functype","mass","charge","atomidx","resid","atomnumber","atomname","atomtype","resname","group"])

# 5. combine info and coord
total=pd.concat([info,coord],axis=1)
total.insert(loc=0,column="index",value=np.arange(1,total.shape[0]+1))
total[["functype","atomidx","resid","atomnumber"]]=total[["functype","atomidx","resid","atomnumber"]].apply(pd.to_numeric)
total[["mass","charge","x","y","z"]]=total[["mass","charge","x","y","z"]].astype(float)

# 6. output
fp=open("output.gro","w")
fp.write("this file was generated by CASEA \n")
fp.write("{} \n".format(total.shape[0]))
for index,row in total.iterrows():
    fp.write(("%5d%-5s%5s%5d%8.3f%8.3f%8.3f \n")%(row["resid"],row["resname"],row["atomname"],row["index"],row["x"],row["y"],row["z"]))
fp.write("{}{}{} \n".format(box[0][0],box[1][1],box[2][2]))
fp.close()

endT=time.time()

print("********************************")
print("tpr2gro is down ")
print("time used is {:3.2f} s ".format(endT-startT))
print("********************************")