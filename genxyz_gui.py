import numpy as np
import pandas as pd
import math
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
 
def select_file1():
    selected_file_path = filedialog.askopenfilename()  # 使用askopenfilename函数选择单个文件
    select_path1.set(selected_file_path) 

def select_file2():
    selected_file_path = filedialog.askopenfilename()  # 使用askopenfilename函数选择单个文件
    select_path2.set(selected_file_path) 
    
def calc():
    try:
        with open(Entry1.get()) as f1:
            m1=f1.readlines()
        xyz1=list(filter(None,[ l.split() for l in m1[2:] if not l.startswith('\n') ]))
        with open(Entry2.get()) as f2:
            m2=f2.readlines()
        xyz2=list(filter(None,[ l.split() for l in m2[2:] if not l.startswith('\n') ]))
        df1=pd.DataFrame(xyz1,columns=["atom","x","y","z"])
        df2=pd.DataFrame(xyz2,columns=["atom","x","y","z"])
        data1=df1[["x","y","z"]].values.astype(np.float64)
        data2=df2[["x","y","z"]].values.astype(np.float64)
        vec=data1[int(Entry3.get())-1]-data2[int(Entry4.get())-1]
        proc_data=data2
        for i in range(1,int(Entry6.get())+1):
            proc_data=np.append(proc_data,data2+vec*float(Entry5.get())*i/math.sqrt(np.square(vec).sum()),axis=0)  
        mm=pd.DataFrame(proc_data,columns=["x","y","z"])
        mm.insert(0,"atom",np.tile(df2["atom"].values,int(Entry6.get())+1))
        # -------------------------------------------------------------------------------- #
        # 输出轨迹
        fp=open("traj.xyz","w")
        for m in range(int(Entry6.get())+1):
            print("{}".format(data1.shape[0]+data2.shape[0]),file=fp)
            print("{}".format(m*float(Entry5.get())),file=fp)
            for i in range(data1.shape[0]):
                print ("{}     {}  {} {}".format(df1.iloc[i]["atom"],df1.iloc[i]["x"],df1.iloc[i]["y"],df1.iloc[i]["z"]),file=fp)
            for j in range(data2.shape[0]):
                print ("{}     {:.8f}  {:.8f} {:.8f}".format(mm.iloc[(m*(data2.shape[0])+j)]["atom"],mm.iloc[(m*(data2.shape[0])+j)]["x"],mm.iloc[(m*(data2.shape[0])+j)]["y"],mm.iloc[(m*(data2.shape[0])+j)]["z"]),file=fp)
        fp.close()
        #-------------------------------------------------------------------------------#
        # 输出gaussian文件
        fg=open("traj.gjf","w")
        for m in range(int(Entry6.get())+1):
            chkname="%chk="+"traj_"+str(m)+".chk"
            print("{}".format(chkname),file=fg)
            print("{}".format("%nproc="+str(Entry9.get())),file=fg)
            print("{}".format("%mem="+str(Entry10.get())),file=fg)
            print("{}".format("#p "+str(Entry11.get())+"/"+str(Entry12.get())+" "+str(Entry14.get())),file=fg)
            print("{}".format("            "),file=fg)
            print("{}".format("traj_"+str(m)+" Structure"),file=fg)
            print("{}".format("            "),file=fg)
            print("{}".format(str(Entry7.get())+" "+str(Entry8.get())),file=fg)
            for i in range(data1.shape[0]):
                print ("{}     {}  {} {}".format(df1.iloc[i]["atom"],df1.iloc[i]["x"],df1.iloc[i]["y"],df1.iloc[i]["z"]),file=fg)
            for j in range(data2.shape[0]):
                print ("{}     {:.8f}  {:.8f} {:.8f}".format(mm.iloc[(m*(data2.shape[0])+j)]["atom"],mm.iloc[(m*(data2.shape[0])+j)]["x"],mm.iloc[(m*(data2.shape[0])+j)]["y"],mm.iloc[(m*(data2.shape[0])+j)]["z"]),file=fg)
            print("{}".format("            "),file=fg)
            # print("{}".format("{}".format("notatoms="+str(Entry13.get()))),file=fg)
            print("{}".format("            "),file=fg)
            if m!=int(Entry6.get()):
                print("{}".format("--link1--"),file=fg)
            else:
                continue
        fg.close()
        messagebox.showinfo("congratulation!","Generating done")      
    except :
        messagebox.showerror("Waring!!!!","please check your input!")

 
if __name__ == '__main__':
    root=Tk()
    root.title("genxyz")
    root.geometry("500x500")
    root.configure(background="#C6DBEF")
    select_path1=StringVar()
    select_path2=StringVar()
    Label1=Label(root,text="xyz1文件路径：").grid(row=0,column=1,padx=10,pady=5)
    Entry1=Entry(root,textvariable=select_path1)
    Entry1.grid(row=0,column=2)
    Button1=Button(root,text="选择xyz文件",command=select_file1).grid(row=0,column=3,padx=10,pady=5)
    Label2=Label(root,text="xyz2文件路径：").grid(row=1,column=1)
    Entry2=Entry(root,textvariable=select_path2)
    Entry2.grid(row=1,column=2)
    Button2=Button(root,text="选择xyz文件",command=select_file2).grid(row=1,column=3)
    atx1=StringVar()
    atx1.set("163")
    Label3=Label(root,text="atomidx1=").grid(row=2,column=1,ipadx=0)
    Entry3=Entry(root,textvariable=atx1)
    Entry3.grid(row=2,column=2,padx=0,pady=5)
    atx2=StringVar()
    atx2.set("1")
    Label4=Label(root,text="atomidx2=").grid(row=2,column=3,padx=0)
    Entry4=Entry(root,textvariable=atx2)
    Entry4.grid(row=2,column=4,padx=0,pady=5)
    stepsize=StringVar()
    stepsize.set("0.05")
    Label5=Label(root,text="stepsize=").grid(row=3,column=1,ipadx=0)
    Entry5=Entry(root,textvariable=stepsize)
    Entry5.grid(row=3,column=2,padx=0,pady=5)
    step=StringVar()
    step.set("50")
    Label6=Label(root,text="step=").grid(row=3,column=3,padx=0)
    Entry6=Entry(root,textvariable=step)
    Entry6.grid(row=3,column=4,padx=0,pady=5)
    charge=StringVar()
    charge.set("0")
    Label7=Label(root,text="charge=").grid(row=4,column=1,ipadx=0)
    Entry7=Entry(root,textvariable=charge)
    Entry7.grid(row=4,column=2,padx=0,pady=5)
    spin=StringVar()
    spin.set("1")
    Label8=Label(root,text="spin=").grid(row=4,column=3,padx=0)
    Entry8=Entry(root,textvariable=spin)
    Entry8.grid(row=4,column=4,padx=0,pady=5)
    nproc=StringVar()
    nproc.set("32")
    Label9=Label(root,text="nproc=").grid(row=5,column=1,ipadx=0)
    Entry9=Entry(root,textvariable=nproc)
    Entry9.grid(row=5,column=2,padx=0,pady=5)
    mem=StringVar()
    mem.set("12GB")
    Label10=Label(root,text="mem=").grid(row=5,column=3,padx=0)
    Entry10=Entry(root,textvariable=mem)
    Entry10.grid(row=5,column=4,padx=0,pady=5)
    method=StringVar()
    method.set("b3lyp")
    Label11=Label(root,text="method=").grid(row=6,column=1,ipadx=0)
    Entry11=Entry(root,textvariable=method)
    Entry11.grid(row=6,column=2,padx=0,pady=5)
    basis=StringVar()
    basis.set("6-31G(d,p)")
    Label12=Label(root,text="basis=").grid(row=6,column=3,padx=0)
    Entry12=Entry(root,textvariable=basis)
    Entry12.grid(row=6,column=4,padx=0,pady=5)
    #fixatom=StringVar()
    #fixatom.set("1,2,16,17,33,34,52,71,78,86,97,111,113,129,130,156,145")
    #Label13=Label(root,text="fixatoms=").grid(row=7,column=1,ipadx=0)
    #Entry13=Entry(root,textvariable=fixatom,width=50)
    #Entry13.grid(row=7,column=2,padx=0,pady=5,columnspan=3)
    other=StringVar()
    other.set("")
    Label14=Label(root,text="other=").grid(row=8,column=1,padx=0)
    Entry14=Entry(root,textvariable=other)
    Entry14.grid(row=8,column=2,padx=0,pady=5)
    
    Button3=Button(text="generator",command=calc)
    Button3.grid(row=9,column=1)
    
    root.mainloop()
    