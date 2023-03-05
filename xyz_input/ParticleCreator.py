from ase import Atoms
from ase.io import write
import numpy as np
one_D_count = 5
# x = np.zeros((one_D_count**2  , 3))
x = np.zeros((one_D_count**2 + (one_D_count+60)*2 +(one_D_count+20)*4-24, 3))
x= x.tolist()
num=0

for i in range(0,one_D_count):
    for j in range(0,one_D_count):        
        x[num][2]=j*2
        x[num][1]=i*2
        num+=1
        
for i in range(-50,one_D_count+10):
    for j in range(0,2):
        x[num][2]=j*2 - 8
        x[num][1]=i*2 
        num+=1
        
for i in range(0,2):
    for j in range(-4,one_D_count+10):
        x[num][2]=j*2 
        x[num][1]=i*2 - 100
        num+=1
        
for i in range(0,2):
    for j in range(-4,one_D_count+10):
        x[num][2]=j*2 
        x[num][1]=i*2 + 120
        num+=1
atom_object = Atoms(positions=x)
write('./xyz_output/particles.xyz', atom_object)
