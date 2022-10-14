import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import sys

file=sys.argv[1]
ss_data=[]
res=[]
with open(file) as f:
    lines=f.readlines()
    lidx=[lines.index(l) for l in lines if l.startswith("  #  RESIDUE AA")]
    data=lines[(lidx[0]+1):]
for i in range(len(data)):   
    if data[i].split()[4] not in ["H","B","E","G","I","P","T","S"]:
        ss_data.append([data[i].split()[3],"T"])
    else:
        ss_data.append([data[i].split()[3],data[i].split()[4]])
    res.append(data[i].split()[3])
ss_dict = {'H': 'H', 'G': 'H', 'I': 'H',
           'B': 'E', 'E': 'E', 'P': 'H',
           'T': 'T', 'S': 'T'}

# 获得不同的二级结构
ss_block=[]
prev_ss=None
for idx,com in enumerate(ss_data):
    reduced_elem = ss_dict.get(com[1], 'T')
    if idx!=0 and idx%50!=0 :
        if reduced_elem != prev_ss:
            ss_block.append([reduced_elem,idx,idx])
            prev_ss = reduced_elem
        ss_block[-1][-1] = idx
    elif idx==0:
        if reduced_elem != prev_ss:
            ss_block.append([reduced_elem,idx,idx])
            prev_ss = reduced_elem
        ss_block[-1][-1] = idx
    else:
        ss_block.append([reduced_elem,idx,(idx%50+1)*50])
        prev_ss = reduced_elem
        ss_block[-1][-1] = idx


reslen=len(ss_data)
fs=9.4
fig=plt.figure(figsize=(10,8),dpi=300,facecolor="white")
ax=fig.add_subplot(111)
# 取消边框
ax.axis("off")
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('top')   
for m,n in enumerate(ss_block):
    ss_type,start,end=n
    if ss_type=="T":
        ax.add_patch(mpl_patches.Rectangle(((start%50)/50-0.001,0.025+(1-((start//50+1))*0.1)),(end-start+1)/50+0.001,0.003,color="gray")) 
        ax.text((start%50)/50+0.0005,0.05+(1-((start//50+1))*0.1)," ".join(res[start:end+1])+" ",size=fs,fontfamily="monospace",fontweight="bold")
        ax.text((start%50)/50,0.072+(1-((start//50+1))*0.1),"|",size=5,color="black")
        ax.text((start%50)/50,0.085+(1-((start//50+1))*0.1),start+1,size=6,color="black")
    elif ss_type=="H":
        n_turns=np.ceil((end-start+1)/1.5)
        x_val=np.linspace((start%50)/50+0.0025,((end+0.95)%50)/50-0.0025,100)
        y_val=(1-((start//50+1))*0.1)+0.1*(-0.4*np.sin(np.linspace(0,n_turns*2*np.pi,100))+1)/4
        ax.plot(x_val,y_val,linewidth=2.4,color="orange",scalex=False, scaley=False)
        ax.text((start%50)/50+0.0005,0.05+(1-((start//50+1))*0.1)," ".join(res[start:end+1])+" ",size=fs,fontfamily="monospace",fontweight="bold")
        ax.text((start%50)/50,0.072+(1-((start//50+1))*0.1),"|",size=5,color="black")
        ax.text((start%50)/50,0.085+(1-((start//50+1))*0.1),start+1,size=6,color="black")
    elif ss_type=="E":
        ax.add_patch(mpl_patches.Rectangle(((start%50)/50,0.01+(1-((start//50+1))*0.1)),(end-start+1)/50-0.0029,0.03,color="lightseagreen"))
        ax.text((start%50)/50+0.0005,0.05+(1-((start//50+1))*0.1)," ".join(res[start:end+1])+" ",size=fs,fontfamily="monospace",fontweight="bold")
        ax.text((start%50)/50,0.072+(1-((start//50+1))*0.1),"|",size=5,color="black")
        ax.text((start%50)/50,0.085+(1-((start//50+1))*0.1),start+1,size=6,color="black")

ax.text(((ss_block[-1][2]+1)%50)/50+0.005,0.019+(1-((ss_block[-1][2]//50+1))*0.1),"C loop",size=fs,fontfamily="monospace",color="gray")
ax.text(((ss_block[-1][2])%50)/50,0.072+(1-((ss_block[-1][2]//50+1))*0.1),"|",size=5,color="black")
ax.text(((ss_block[-1][2])%50)/50,0.085+(1-((ss_block[-1][2]//50+1))*0.1),reslen,size=6,color="black")

plt.savefig("ss.png",dpi=300)

