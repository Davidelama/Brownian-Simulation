import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
from mpl_toolkits.mplot3d import Axes3D
import sys, getopt

def main(step=99):
  par=np.genfromtxt("positions_step"+str(step)+".dat",max_rows=1)
  L=par[0]
  radius=par[1]/2
  rc=par[2]
  data=np.genfromtxt("positions_step"+str(step)+".dat",skip_header=1)
  x=data[:,0]
  y=data[:,1]
  z=data[:,2]
  
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(x, y, z, c='r', marker='o')
  ax.set_xlim(-L/2,L/2)
  ax.set_ylim(-L/2,L/2)
  ax.set_zlim(-L/2,L/2)
  plt.show()
  fig.savefig(r'Positions'+str(step)+'.pdf', bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(step=int(sys.argv[1]))
    else:
        sys.exit('usage: pos_plot.py <step>')