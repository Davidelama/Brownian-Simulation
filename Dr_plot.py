import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import pandas as pd


par=np.genfromtxt("Dr.dat",max_rows=1)
data=np.genfromtxt("Dr.dat",skip_header=1)#, skip_footer=1
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



fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(14,6))
for i in range(3):
  ax[0].errorbar(time, corr[i,:],yerr=dcorr[i,:]) #a BEAN
#ax.set_ylim(0,)
ax[0].set_yscale("log")
ax[0].set_xscale("log")
ax[0].set_title("Mean square distance as function of time")
ax[0].set_xlabel("t")
ax[0].set_ylabel(r"$\langle \Delta r^2\rangle$")

D=np.zeros((3,Ncorr-1))
dD=np.zeros((3,Ncorr-1))

htime=np.zeros((Ncorr-1))


for i in range(3):
  for j in range(Ncorr-1):
    htime[j]=(time[j+1]+time[j])*.5
    D[i,j]=(corr[i,j+1]-corr[i,j])/(time[j+1]-time[j])*.5
    dD[i,j]=(np.sqrt(dcorr[i,j+1]**2+dcorr[i,j]**2))/(time[j+1]-time[j])*.5
  ax[1].errorbar(htime, D[i,:],yerr=dD[i,:]) #a BEAN
  
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_title("Dispersion as function of time")
ax[1].set_xlabel("t")
ax[1].set_ylabel(r"D")  

plt.show()
fig.savefig(r'Dr.pdf', bbox_inches='tight')