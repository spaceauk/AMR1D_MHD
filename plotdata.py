import numpy as np
import matplotlib.pyplot as plt

# Read data from the file
data = np.loadtxt('data/finalresult.dat') 

# Sort the data based on the first column
sorted_data = data[data[:,0].argsort()]

# Extract x values and prim variables
x_values = sorted_data[:,0]
prim_variables = sorted_data[:, 1:]
num_prim_variables = prim_variables.shape[1]
print("Number of variables: "+str(num_prim_variables))

# Define markers for different prim variables
if num_prim_variables==3:
    markers = ['.', '.', '.']
    ylabel=['Density','Velocity','Pressure']
    vert_size=4
else:
    markers = ['.', '.', '.','.','.','.']
    ylabel=['Density','Vx','Pressure','Vy','Bx','By']
    vert_size=8

# Plot x with respect to each prim variable
num_prim_variables = prim_variables.shape[1]
plt.figure(figsize=(14,vert_size))
for i in range(num_prim_variables):
    if num_prim_variables==3:
        plt.subplot(1,3,i+1)
    else:
        plt.subplot(2,3,i+1)
    plt.plot(x_values, prim_variables[:, i], marker=markers[i],label=f'prim_{i+1}')
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel(ylabel[i])
    if i==1:
        plt.title('Plot of Prim Variables vs. x')
print("Mesh size (Python): "+str(len(x_values)))
plt.show()

