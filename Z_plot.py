import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import pandas as pd


par=np.genfromtxt("Z_corr.dat",max_rows=1)
data=np.genfromtxt("Z_corr.dat",skip_header=1, skip_footer=1)
Nblock=int(par[0])
Nlevel=int(par[1])
dt=par[2]
Nmeas=len(data[:,0])
Ncorr=(Nblock-1)*Nlevel
raw_data={"T":data[:,0],"xmeas":data[:,1], "ymeas":data[:,2], "zmeas":data[:,3]}
df=pd.DataFrame(raw_data,columns=["T","xmeas","ymeas","zmeas"])

corr=np.zeros((3,Ncorr))
dcorr=np.zeros((3,Ncorr))
time=np.zeros((Ncorr))

for i in range(Ncorr):
  time[i]=(i+1)%(Nblock-1)
  if time[i]==0:
    time[i]=(Nblock-1)/Nblock
  time[i]*=dt*Nblock**((i+1)//(Nblock-1))
  cond=df["T"]==i
  obj=df[cond].to_numpy()
  #print(i, len(obj))
  for j in range(3):
    corr[j,i]=np.mean(obj[:,j+1])
    dcorr[j,i]=np.std(obj[:,j+1])/np.sqrt(len(obj[:,j+1])) # we are "forgetting" autocorrelation



fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,6))
for i in range(3):
  ax.errorbar(time, corr[i,:],yerr=dcorr[i,:]) #a BEAN
ax.set_ylim(0,)
ax.set_xscale("log")
ax.set_title("Velocity correlation as function of time")
ax.set_xlabel("t")
ax.set_ylabel(r"$\langle v(t)v(0)\rangle$")
plt.show()
fig.savefig(r'Z_total.pdf', bbox_inches='tight')