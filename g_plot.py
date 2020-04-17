import numpy as np
import matplotlib.pyplot as plt
import sys, getopt

def main(step=99):
  
  par=np.genfromtxt("positions_step"+str(step)+".dat",max_rows=1)
  
  data=np.genfromtxt("g_step"+str(step)+".dat",skip_header=1)
  fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8,6))
  ax.step(data[:,0]+data[1,0], data[:,1]) #a BEAN
  ax.set_xlim(0,)
  ax.set_ylim(0,)
  plt.show()
  fig.savefig(r'g'+str(step)+'.pdf', bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(step=int(sys.argv[1]))
    else:
        sys.exit('usage: g_plot.py <step>')