import pandas as pd
import numpy as np
import math
from MDAnalysis.coordinates.XTC import XTCReader
import matplotlib.pyplot as plt
import matplotlib
import sys

class Traj:
    "Reads trajectory from .xtc file"
                
    def __init__(self,fn):
        '''
        # Paramters \n  
        fn: traj filename
        b: start time ps (default -1)
        e: end time ps(default -1)
        dt: skip time ps (default 2 only even numbers)
        nc: cores of cpu to use (default 4)
        '''
        # dt即时间间隔只能为偶数
        self.fn=fn
        self.traj=self.readTraj()
        self.coord=self.getCoord()

        
    def readTraj(self):
        return XTCReader(self.fn)
        
    # def getBox(self):
    #     return [self.traj.box[i][i]  for i in range(3)]
    
    def getAtomNum(self):
        return self.traj.n_atoms
    
    def getCoord(self):
        # 0对应于vmd中的第一帧
        atomCoord={}
        x=np.array([[0 for i in range(1,4)] for j in range(1,self.getAtomNum()+1)])
        for t in self.traj.trajectory:
            x=np.array(t)
            atomCoord[str(t.time)]=x
        n_series=pd.Series(atomCoord)
        tmp=pd.concat([(pd.DataFrame(n_series[i],columns=[[str(n_series.keys()[i]),str(n_series.keys()[i]),str(n_series.keys()[i])],["x","y","z"]])) for i in range(n_series.size)],axis=1)
        print("input file contains {} frames  from {} ps to {} ps dt {}ps".format(int(tmp.shape[1]/3),n_series.keys()[0],n_series.keys()[-1],float(n_series.keys()[1])-float(n_series.keys()[0])))
        # print(tmp)
        return tmp
    
class gro:
    "读取gro文件"
    
    def __init__(self,fn):
        self.fn=fn 
        self.atomInfo=self.readGro()
        
    def readGro(self):
        gro=open(self.fn,"r")
        atomSet=[]
        for line in gro.readlines():
            tmp=[]
            if len(line.split())>3:
                # 1 resnum 2 resn 3 aton name 4 atom num 5x 6y 7z 8vx 9vy 10vz
                tmp.append(line[0:5].replace(" ",""))
                tmp.append(line[5:10].replace(" ",""))
                tmp.append(line[10:15].replace(" ",""))
                tmp.append(line[15:20].replace(" ",""))
                tmp.append(line[20:28].replace(" ",""))
                tmp.append(line[28:36].replace(" ",""))
                tmp.append(line[36:44].replace(" ",""))
                tmp.append(line[44:52].replace(" ",""))
                tmp.append(line[52:60].replace(" ",""))
                tmp.append(line[60:68].replace(" ",""))
            if tmp!=[]:
                atomSet.append(tmp)
        return pd.DataFrame(atomSet,columns=["resNum","resName","atomName","atomNum","X","Y","Z","Vx","Vy","Vz"])
    
class combine():
    """
    读取gro文件和xtc文件并整合
    """
    def __init__(self,grofile,xtcfile):
        '''
        # Paramters \n  
        fn: traj filename
        b: start frame (default -1)
        e: end frame (default -1)
        dt: skip frames (default 2 only even numbers)
        nc: cores of cpu to use (default 4)
        '''
        self.grofile=grofile
        self.xtcfile=xtcfile
        
    def gro_xtc(self):
        groFile=gro(self.grofile).readGro()
        xtcFile=Traj(self.xtcfile).coord
        xtcFile.insert(0,"atomNum",groFile["atomNum"])
        xtcFile.insert(1,"atomName",groFile["atomName"])
        xtcFile.insert(2,"resNum",groFile["resNum"])
        xtcFile.insert(3,"resName",groFile["resName"])
        return xtcFile   
    
def resi_sele(sheet,sele):
    """
    #Paramters
    sele example: `resi 1:5`
                  `resi 1:5,10`
                  `resi 1,2,3`
    """
    sele_name=list(filter(None,sele.split(" ")))
    if sele_name[0].lower()=="resi":
        tmp=[]
        for i in sele_name[1].split(","):
            if len(i.split(":"))==2:
                tmp.append(list(map(str,range(int(i.split(":")[0]),int(i.split(":")[1])+1))))
            else:
                tmp.append([i])
            return sheet[sheet["resNum"].isin(sum(tmp,[]))]
    else:
        print("please enter the correct format(resi 1:5 or resi 1,2) ")
        

def resn_sele(sheet,sele):
    """
    #Paramters
    only return existed resname
    sele example: `resn xxx,xxx`
    """
    sele_name=list(filter(None,sele.split(" ")))
    if sele_name[0].lower()=="resn":
        tmp=[]
        for i in sele_name[1].split(","):
            tmp.append([i]) 
        return sheet[sheet["resName"].isin(sum(tmp,[]))] 
    else:
        print("please enter the correct format(resn xxx) ") 
    

def atom_sele(sheet,sele):
    """
    #Paramters
    only return existed atomname
    sele example: `name xxx,xxx`
    if use `name C*` return all C atom
    """
    sele_name=list(filter(None,sele.split(" ")))
    tmp=[]
    if sele_name[0].lower()=="name":
        for i in sele_name[1].split(","):
            if i.endswith("*"):
                tmp.append(sheet[sheet["atomName"].str.contains(i[0])]["atomName"].unique().tolist())   
            else:
                tmp.append([i])
        return sheet[sheet["atomName"].isin(sum(tmp,[]))]  
    else:
        print("please enter the correct format(name C1,C2 or name C*) ")


m=combine(sys.argv[1],sys.argv[2]).gro_xtc()

CA=atom_sele(m,"name CA").reset_index()

tmp=[]
for i in range(int((CA.shape[1]-5)/3)):
    tmp.append(np.array(CA.iloc[:,(5+i*3):(8+i*3)]))
allxyz=np.array(tmp)
# tmp储存格式帧，原子，x，y，z

# kabschalgorithm 点云算法修正坐标 fit
# fitxyz
for n in range(allxyz.shape[0]):
    cercoord_A=np.mean(np.matrix(allxyz[0]),axis=0)
    cercoord_B=np.mean(np.matrix(allxyz[n]),axis=0)
    A_c=allxyz[0]-cercoord_A
    B_c=allxyz[n]-cercoord_B
    H=B_c.T*A_c
    U,S,VT=np.linalg.svd(H)
    d=np.sign(np.linalg.det(VT.T*U.T))
    diag=np.diag([1,1,d])
    R=VT.T * diag *  U.T
    # 坐标以A为主
    T=-R*cercoord_B.T+cercoord_A.T
    # 修正之后的B坐标
    B_calc=(R*allxyz[n].T+T).T
    allxyz[n]=B_calc

atomcoord_mean=[]
for i in range(allxyz.shape[1]):
    x=[]
    y=[]
    z=[]
    for j in range(allxyz.shape[0]):
        x.append(allxyz[j][i][0])
        y.append(allxyz[j][i][1])
        z.append(allxyz[j][i][2])
    atomcoord_mean.append([np.mean(x),np.mean(y),np.mean(z)])
delta0=np.array(atomcoord_mean)   

## new method to calculate the covariance and correlation matrix
## Two methods produce numbers with deviation less than 0.00001
offset = allxyz
# print(offset.shape) # (10001, 130, 3)
# np.cov for x, y, and z, then sum
## 这里对数据重新处理了一下，把坐标的三维拆开，然后重新处理成(原子，时间帧)的形式
## 例如offset_X是包含了130个列表的变量，每一个列表里面保存了一个原子的X坐标随时间变化的10001个数据
offset_X, offset_Y, offset_Z = offset[:, :, 0].T, offset[:, :, 1].T, offset[:, :, 2].T
# print(offset_X.shape) # (130, 10001)
## 对每一个维度的数据做协方差矩阵，得到的是一个（130，130）的矩阵
## 也即每一个维度上，原子之间的协方差矩阵
# np.cov for x, y, and z, then sum
covariance_X = np.cov(offset_X, ddof=0) # /n not /(n‐1)
covariance_Y = np.cov(offset_Y, ddof=0)
covariance_Z = np.cov(offset_Z, ddof=0)
## 将三个维度相加，得到最终的协方差矩阵
covariance = covariance_X + covariance_Y + covariance_Z
# print(covariance_Z.shape) # （130， 130）
# print(covariance.shape) # （130， 130）
## 将协方差矩阵转换为互相关矩阵
corr=np.zeros((atom_number,atom_number))
for i in range(0,atom_number):
    for j in range(0,atom_number):
        corr[i,j] = covariance[i,j]/np.sqrt(covariance[i,i]*covariance[j,j])

# 计算差值         
# 以平均原子坐标为基准
# np.set_printoptions(precision=5)
# delta=allxyz-delta0
# delta_r=np.empty((delta.shape[1],delta.shape[1],delta.shape[0],))
# for n in range(delta.shape[0]):
#     for i in range(delta.shape[1]):
#         for j in range(delta.shape[1]):
#             delta_r[i,j,n]=np.dot(delta[n][i][:],delta[n][j][:])
# cij=delta_r           

# # 先时间平均
# deltai=np.empty((delta.shape[1]))
# for i in range(delta.shape[1]):
#     deltai[i]=np.sqrt(np.mean(cij[i,i,:]))
        
# Cij=np.empty((delta.shape[1],delta.shape[1]))
# for i in range(delta.shape[1]):
#     for j in range(delta.shape[1]):
#         Cij[i,j]=np.mean(cij[i,j,:])/(deltai[i]*deltai[j])
        
Cout=corr


# 使用contour+cp=ontour绘图
# 使用contour+cp=ontour绘图
plt.figure(figsize=(14,12))
cmap1 = matplotlib.colors.ListedColormap(['#ff80ff',"#ffabff","#ffd4ff","#fffcff","#d4ffff","#a8ffff","#80ffff"])
norm=matplotlib.colors.BoundaryNorm([-1.0,-0.75,-0.5,-0.25,0.25,0.5,0.75,1.0],cmap1.N)
x=np.array(range(Cij.shape[0]))
y=np.array(range(Cij.shape[1]))
z=Cij
plt.rcParams['font.size'] = 16
plt.contourf(x,y,z,cmap=cmap1,levels=[-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1],norm=norm)
plt.colorbar()
plt.contour(x,y,z,colors="grey",levels=[-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1],linewidths=0.8,linestyles="solid")
ax=plt.gca()
font1 = {'family' : 'Arial',
         'weight' : 'normal',
         'size'   : 24,
         }
ax.set_xlabel("Residue No.",font1,labelpad=15.0)
ax.set_ylabel("Residue No.",font1,labelpad=15.0),
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(25))
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(25))
ax.set_title("Residue Cross Correlation",fontweight="bold",fontsize=26,pad=20)
plt.savefig("DCCM.png")
