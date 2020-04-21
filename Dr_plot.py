import matplotlib as mpl
#mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import pandas as pd


par=np.genfromtxt("Dr.dat",max_rows=1)
data=np.genfromtxt("Dr.dat",skip_header=1, skip_footer=1)#
Nblock=int(par[0])
Nlevel=int(par[1])
dt=par[2]
Nmeas=len(data[:,0])
Ncorr=(Nblock-1)*Nlevel
raw_data={"T":data[:,0],"xmeas":data[:,1], "ymeas":data[:,2], "zmeas":data[:,3]}
df=pd.DataFrame(raw_data,columns=["T","xmeas","ymeas","zmeas"])

corr=np.zeros((Ncorr))
corr2=np.zeros((Ncorr))
alpha=np.zeros((Ncorr))
dcorr=np.zeros((Ncorr))
time=np.zeros((Ncorr))

for i in range(Ncorr):
  time[i]=(i+1)%(Nblock-1)
  if time[i]==0:
    time[i]=(Nblock-1)/Nblock
  time[i]*=dt*Nblock**((i+1)//(Nblock-1))
  cond=df["T"]==i
  obj=df[cond].to_numpy()
  #print(i, len(obj))
  temp=np.concatenate((obj[:,1],obj[:,2],obj[:,3]))
  corr[i]=np.mean(temp)
  dcorr[i]=np.std(temp)/np.sqrt(len(temp)) # we are "forgetting" autocorrelation
  corr2[i]=np.mean(temp*temp)
  #dcorr2[i]=np.std(temp*temp)/np.sqrt(len(temp))
  alpha[i]=1-corr2[i]/(3*corr[i]*corr[i])



fig, ax = plt.subplots(nrows=1,ncols=3,figsize=(14,6))
ax[0].errorbar(time, corr,yerr=dcorr) #a BEAN
#ax.set_ylim(0,)
ax[0].set_yscale("log")
ax[0].set_xscale("log")
ax[0].set_title("Mean square distance as function of time")
ax[0].set_xlabel("t")
ax[0].set_ylabel(r"$\langle \Delta r^2\rangle$")

D=np.zeros((Ncorr-1))
dD=np.zeros((Ncorr-1))

htime=np.zeros((Ncorr-1))


for j in range(Ncorr-1):
  htime[j]=(time[j+1]+time[j])*.5
  D[j]=(corr[j+1]-corr[j])/(time[j+1]-time[j])*.5
  dD[j]=(np.sqrt(dcorr[j+1]**2+dcorr[j]**2))/(time[j+1]-time[j])*.5
ax[1].errorbar(htime, D,yerr=dD) #a BEAN
  
ax[1].set_yscale("log")
ax[1].set_xscale("log")
ax[1].set_title("Dispersion as function of time")
ax[1].set_xlabel("t")
ax[1].set_ylabel(r"D")  

ax[2].plot(time, alpha) #a BEAN
  
ax[2].set_xlabel("t")
ax[2].set_ylabel(r"$\alpha$") 
ax[2].set_xscale("log")

plt.show()
fig.savefig(r'Dr.pdf', bbox_inches='tight')
