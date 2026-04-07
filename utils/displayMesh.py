# Use grid.dat file as input argument
#
# $ python3 displayMesh.py relative/path/to/grid.dat
import matplotlib.pyplot as pl
import numpy as np
import sys

cl = 'white'

filename = sys.argv[1]
fc=open(filename,'r')
line=fc.readline()

nodes= int(line.split()[0])
neles = int(line.split()[1])

A = 1
B = neles + 1 - nodes
C = neles
mn = np.roots([A,B,C])

x=[]
y=[]
for i in range(nodes):
 line=fc.readline()
 x.append(float(line.split()[1]))
 y.append(float(line.split()[2]))

pl.style.use('dark_background')
pl.scatter(x, y, color=cl, marker='.', s=0.8)
pl.axis('equal')
pl.axis([min(x), max(x)*0.9, 0.0, max(y)*1.1])
pl.xlabel('x (mm)')
pl.ylabel('y (mm)')
pl.title(filename+' points: '+str(int(mn[0]+1))+'x'+str(int(mn[1]+1))+' grid')

for i in range(neles):
 line = fc.readline()
 n1 = int(line.split()[1])-1
 n2 = int(line.split()[2])-1
 n3 = int(line.split()[3])-1
 n4 = int(line.split()[4])-1
 pl.plot([x[n1],x[n2],x[n3],x[n4],x[n1]], [y[n1],y[n2],y[n3],y[n4],y[n1]], color = cl, linewidth=0.7)
 
fc.close()
pl.show()



