import numpy as np
import pandas as pd
import math


# 初始参数设置
fn1="1.xyz"
fn2="2.xyz"
atomidx1=163
atomidx2=1
stepsize=0.05
step=50
charge=0
spin=1
nproc=32
mem="12GB"
method="b3lyp"
basis="6-31g(d,p)"
#fixatoms="1,2,16,17,33,34,52,71,78,86,97,111,113,129,130,156,145"
other=""
#------------------------------------------------------------------------------#
with open(fn1) as f1:
    m1=f1.readlines()
xyz1=list(filter(None,[ l.split() for l in m1[2:] if not l.startswith('\n') ]))
with open(fn2) as f2:
    m2=f2.readlines()
xyz2=list(filter(None,[ l.split() for l in m2[2:] if not l.startswith('\n') ]))
df1=pd.DataFrame(xyz1,columns=["atom","x","y","z"])
df2=pd.DataFrame(xyz2,columns=["atom","x","y","z"])
data1=df1[["x","y","z"]].values.astype(np.float64)
data2=df2[["x","y","z"]].values.astype(np.float64)
vec=data1[atomidx1-1]-data2[atomidx2-1]
proc_data=data2
for i in range(1,step+1):
    proc_data=np.append(proc_data,data2+vec*stepsize*i/math.sqrt(np.square(vec).sum()),axis=0)  
mm=pd.DataFrame(proc_data,columns=["x","y","z"])
mm.insert(0,"atom",np.tile(df2["atom"].values,step+1))
# -------------------------------------------------------------------------------- #
# 输出轨迹
fp=open("traj.xyz","w")
for m in range(step+1):
    print("{}".format(data1.shape[0]+data2.shape[0]),file=fp)
    print("{}".format(m*stepsize),file=fp)
    for i in range(data1.shape[0]):
        print ("{}     {}  {} {}".format(df1.iloc[i]["atom"],df1.iloc[i]["x"],df1.iloc[i]["y"],df1.iloc[i]["z"]),file=fp)
    for j in range(data2.shape[0]):
        print ("{}     {:.8f}  {:.8f} {:.8f}".format(mm.iloc[(m*(data2.shape[0])+j)]["atom"],mm.iloc[(m*(data2.shape[0])+j)]["x"],mm.iloc[(m*(data2.shape[0])+j)]["y"],mm.iloc[(m*(data2.shape[0])+j)]["z"]),file=fp)
fp.close()
#-------------------------------------------------------------------------------#
# 输出gaussian文件
fg=open("traj.gjf","w")
for m in range(step+1):
    chkname="%chk="+"traj_"+str(m)+".chk"
    print("{}".format(chkname),file=fg)
    print("{}".format("%nproc="+str(nproc)),file=fg)
    print("{}".format("%mem="+str(mem)),file=fg)
    print("{}".format("#p "+method+"/"+basis+" "+other),file=fg)
    print("{}".format("            "),file=fg)
    print("{}".format("traj_"+str(m)+" Structure"),file=fg)
    print("{}".format("            "),file=fg)
    print("{}".format(str(charge)+" "+str(spin)),file=fg)
    for i in range(data1.shape[0]):
        print ("{}     {}  {} {}".format(df1.iloc[i]["atom"],df1.iloc[i]["x"],df1.iloc[i]["y"],df1.iloc[i]["z"]),file=fg)
    for j in range(data2.shape[0]):
        print ("{}     {:.8f}  {:.8f} {:.8f}".format(mm.iloc[(m*(data2.shape[0])+j)]["atom"],mm.iloc[(m*(data2.shape[0])+j)]["x"],mm.iloc[(m*(data2.shape[0])+j)]["y"],mm.iloc[(m*(data2.shape[0])+j)]["z"]),file=fg)
    print("{}".format("            "),file=fg)
    #print("{}".format("{}".format("notatoms="+fixatoms)),file=fg)
    print("{}".format("            "),file=fg)
    if m!=step:
        print("{}".format("--link1--"),file=fg)
    else:
        continue
fg.close()