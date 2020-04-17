import numpy as np
import matplotlib.pyplot as plt
import sys, getopt

par=np.genfromtxt("sf_total.dat",max_rows=1)

data=np.genfromtxt("sf_total.dat",skip_header=1)
Tot=len(data[0,:])
Nmeas=len(data[:,0])
q=np.zeros((Tot))
g_avgs=np.zeros((Tot,2))
for i in range(Tot):
  g_avgs[i,0]=np.mean(data[:,i])
  g_avgs[i,1]=np.std(data[:,i])/np.sqrt(Nmeas)

for i in range(Tot):
  q[i]=(i+1)*par

#print(g_avgs)
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,6))
ax.errorbar(q, g_avgs[:,0],yerr=g_avgs[:,1],linestyle="") #a BEAN
ax.set_xlabel("q")
ax.set_ylabel("S")
ax.set_title("Structure factor as function of momenta")
plt.show()
fig.savefig(r'sf.pdf', bbox_inches='tight')