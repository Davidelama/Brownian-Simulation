import numpy as np
import matplotlib.pyplot as plt
import sys, getopt

par=np.genfromtxt("g_total.dat",max_rows=1)

data=np.genfromtxt("g_total.dat",skip_header=1)
Nhist=len(data[0,:])
Nmeas=len(data[:,0])
r=np.zeros((Nhist))
g_avgs=np.zeros((Nhist,2))
for i in range(Nhist):
  g_avgs[i,0]=np.mean(data[:,i])
  g_avgs[i,1]=np.std(data[:,i])#/np.srrt(Nmeas)

for i in range(Nhist):
  r[i]=(i+.5)*par
  
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,6))
ax.errorbar(r, g_avgs[:,0],yerr=g_avgs[:,1],xerr=par*.5*np.ones(Nhist),linestyle="") #a BEAN
ax.set_xlim(0,)
ax.set_ylim(0,)
plt.show()
fig.savefig(r'g_total.pdf', bbox_inches='tight')