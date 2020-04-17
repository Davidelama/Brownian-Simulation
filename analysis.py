import numpy as np
import matplotlib.pyplot as plt
import sys, getopt

def main():
  
  par=np.genfromtxt("info.dat", skip_header=1)
  L=par[0]
  N=par[1]
  T=par[2]
  dt=par[3]
  sigma=par[4]
  rho=par[5]
  
  
  data=np.genfromtxt("output.dat")
  fig1, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(8,6))
  Ek=data[:,0]
  Ek2=Ek*Ek
  U=data[:,1]
  E=Ek+U
  Cv=3*N/2*(1-2/(3*N*T*T)*(np.mean(Ek2)-np.mean(Ek)**2))
  print("Cv={0:.2f}".format(Cv))
  ax1.plot(range(len(Ek)),Ek,label="$E_k={0:.2f} \pm {1:.2f}$".format(np.mean(Ek),np.std(Ek)))
  ax1.plot(range(len(U)),U,label="$U={0:.2f} \pm {1:.2f}$".format(np.mean(U),np.std(U)))
  ax1.plot(range(len(E)),E,label="$E={0:.4f} \pm {1:.4f}$".format(np.mean(E),np.std(E)))
  ax1.set_xlim(0,)
  ax1.set_ylim(0,)
  ax1.legend()
  fig1.savefig(r'Energy.pdf', bbox_inches='tight')
  
  fig2, ax2 = plt.subplots(nrows=1,ncols=1,figsize=(8,6))
  P=data[:,2]
  ax2.plot(range(len(P)),P,label="$P={0:.2f} \pm {1:.2f}$".format(np.mean(P),np.std(P)))
  ax2.set_xlim(0,)
  ax2.set_ylim(0,)
  ax2.legend()
  fig2.savefig(r'Pressure.pdf', bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main();#main(step=int(sys.argv[1]))
    else:
        sys.exit('usage: g_plot.py <step>')