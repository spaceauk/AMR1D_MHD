import numpy as np
import matplotlib.pyplot as plt

# Read data from the file
cdata = np.loadtxt('data/char_curves.dat')

# Re-arrange the data into matrix for contour plotting
numx=0
while cdata[numx,0]==0:
    numx+=1
lent=len(cdata[:,0])
numt=int(lent/numx)
print(f"Number of points along x={numx} and t={numt}")
Jinv=np.zeros((numt,numx,4)) # 0=x-loc, 1=Jplus, 2=Jminus, 3=J0=entropy
Jinv_sort=np.zeros((numt,numx,3))
tdomain=np.zeros((numt,numx))
for t in range(lent):
    Jinv[int(t/numx),int(t%numx),0:3]=cdata[t,1:4]    
    tdomain[int(t/numx),:]=cdata[t,0]
for t in range(numt):
    Jinv[t] = Jinv[t, Jinv[t, :, 0].argsort()]


# Plot x with respect to each prim variable
vert_size=4
ylabel=['J+','J-','s']
plt.figure(figsize=(14,vert_size))
for i in range(3):
    plt.subplot(1,3,i+1)
    plt.contourf(Jinv[:,:,0],tdomain,Jinv[:,:,1+i],levels=50)
    plt.xlabel('x')
    plt.ylabel(ylabel[i])
    if i==1:
        plt.title('Contour plot of Riemann invariants (x vs t)')
plt.show()

