import optparse
import numpy as np
import pandas as pd
import math

usage="Usage: %prog [options] arg1 arg2..."
parser=optparse.OptionParser(usage,version="%prog 1.0")
parser.add_option("-1","--filename1",dest="fn1",default="1.xyz",
                  help="one of the file")
parser.add_option("-2","--filename2",dest="fn2",default="2.xyz",
                  help="another file")
parser.add_option("-3","--atomindex1",dest="atomidx1",default=163,
                  help="atom1 index ")
parser.add_option("-4","--atomindex2",dest="atomidx2",default=1,
                  help="atom2 index ")
parser.add_option("-d","--stepsize",dest="stepsize",default=0.05,
                  help="stepsize ")
parser.add_option("-s","--step",dest="step",default=50,
                  help="step")
parser.add_option("-c","--charge",dest="charge",default=0,
                  help="charge ")
parser.add_option("-i","--spin",dest="spin",default=1,
                  help="spin")
parser.add_option("-n","--nproc",dest="nproc",default=32,
                  help="nproc")
parser.add_option("-m","--memory",dest="mem",default="12GB",
                  help="memory")
parser.add_option("-e","--method",dest="method",default="b3lyp",
                  help="method")
parser.add_option("-b","--basis",dest="basis",default="6-31G(d,p)",
                  help="basis-set")
#parser.add_option("-f","--fixatoms",dest="fixatoms",default="1,2,16,17,33,34,52,71,78,86,97,111,113,129,130,156,145",
#                  help="fixatoms")
parser.add_option("-o","--other",dest="other",default="",
                  help="other command lines")
options,args=parser.parse_args()

# fn1="1.xyz"
# fn2="2.xyz"
# atomidx1=163
# atomidx2=1
# stepsize=0.05
# step=50
# charge=0
# spin=1
# nproc=32
# mem="12GB"
# method="b3lyp"
# basis="6-31g(d,p)"
# fixatoms="1,2,16,17,33,34,52,71,78,86,97,111,113,129,130,156,145"
# other=""

with open(options.fn1) as f1:
    m1=f1.readlines()
xyz1=list(filter(None,[ l.split() for l in m1[2:] if not l.startswith('\n') ]))

with open(options.fn2) as f2:
    m2=f2.readlines()
xyz2=list(filter(None,[ l.split() for l in m2[2:] if not l.startswith('\n') ]))

df1=pd.DataFrame(xyz1,columns=["atom","x","y","z"])
df2=pd.DataFrame(xyz2,columns=["atom","x","y","z"])
data1=df1[["x","y","z"]].values.astype(np.float64)
data2=df2[["x","y","z"]].values.astype(np.float64)
vec=data1[int(options.atomidx1)-1]-data2[int(options.atomidx2)-1]
proc_data=data2
for i in range(1,int(options.step)+1):
    proc_data=np.append(proc_data,data2+vec*float(options.stepsize)*i/math.sqrt(np.square(vec).sum()),axis=0)
mm=pd.DataFrame(proc_data,columns=["x","y","z"])
mm.insert(0,"atom",np.tile(df2["atom"].values,int(options.step)+1))
# 输出轨迹
fp=open("traj.xyz","w")
for m in range(int(options.step)+1):
    print("{}".format(data1.shape[0]+data2.shape[0]),file=fp)
    print("{}".format(m*float(options.stepsize)),file=fp)
    for i in range(data1.shape[0]):
        print ("{}     {}  {} {}".format(df1.iloc[i]["atom"],df1.iloc[i]["x"],df1.iloc[i]["y"],df1.iloc[i]["z"]),file=fp)
    for j in range(data2.shape[0]):
        print ("{}     {:.8f}  {:.8f} {:.8f}".format(mm.iloc[(m*(data2.shape[0])+j)]["atom"],mm.iloc[(m*(data2.shape[0])+j)]["x"],mm.iloc[(m*(data2.shape[0])+j)]["y"],mm.iloc[(m*(data2.shape[0])+j)]["z"]),file=fp)
fp.close()

# 输出gaussian文件
fg=open("traj.gjf","w")
for m in range(int(options.step)+1):
    chkname="%chk="+"traj_"+str(m)+".chk"
    print("{}".format(chkname),file=fg)
    print("{}".format("%nproc="+str(options.nproc)),file=fg)
    print("{}".format("%mem="+str(options.mem)),file=fg)
    print("{}".format("#p  "+options.method+"/"+options.basis+options.other),file=fg)
    print("{}".format("            "),file=fg)
    print("{}".format("traj_"+str(m)+" Structure"),file=fg)
    print("{}".format("            "),file=fg)
    print("{}".format(str(options.charge)+" "+str(options.spin)),file=fg)
    for i in range(data1.shape[0]):
        print ("{}     {}  {} {}".format(df1.iloc[i]["atom"],df1.iloc[i]["x"],df1.iloc[i]["y"],df1.iloc[i]["z"]),file=fg)
    for j in range(data2.shape[0]):
        print ("{}     {:.8f}  {:.8f} {:.8f}".format(mm.iloc[(m*(data2.shape[0])+j)]["atom"],mm.iloc[(m*(data2.shape[0])+j)]["x"],mm.iloc[(m*(data2.shape[0])+j)]["y"],mm.iloc[(m*(data2.shape[0])+j)]["z"]),file=fg)
    print("{}".format("            "),file=fg)
    #print("{}".format("{}".format("notatoms="+options.fixatoms)),file=fg)
    print("{}".format("            "),file=fg)
    if m!=int(options.step):
        print("{}".format("--link1--"),file=fg)
    else:
        continue
fg.close()