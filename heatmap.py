#!/home/florian/anaconda3/bin/python
import matplotlib.pyplot as plt
import numpy as np


i1, i2, force, energycl, energyq = np.loadtxt('principale_component_e0_e1_smooth.tex', unpack= True)

e0 = np.loadtxt('e0_list.txt')
e1 = np.loadtxt('e1_list.txt')
E0 = np.loadtxt('e0_list_smooth.txt')
E1 = np.loadtxt('e1_list_smooth.txt')

energynm = energyq - energycl

#Definition of length of each conformation of PCA

ei = e0
ej = e1

c = len(ei)// 69
        

l = len(ei) - 60*c
d = len(ei) - 50*c
e = len(ei) - 40*c
f = len(ei) - 30*c
p = len(ei) - 20*c

t1_list = []
t2_list = []

for c in range(e, f):
    t1_list.append(ei[c])
    t2_list.append(ej[c])
    
for cc in range(p , len(ei)):
    t1_list.append(ei[cc])
    t2_list.append(ej[cc])
    

x_1, y_1 = ei[0:l], ej[0:l]
x_2, y_2 = ei[l:d], ej[l:d]
x_3, y_3 = ei[d:e], ej[d:e]
x_4, y_4 = t1_list, t2_list
x_5, y_5 = ei[f:p], ej[f:p]


plt.figure()
lw = 0.5
ps = 8.
plt.scatter(x_1, y_1, c='green', s=ps, marker='o', label='1', linewidths=lw, edgecolors="w")
plt.scatter(x_2, y_2, c='purple', s=ps, marker='o', label='2', linewidths=lw, edgecolors="w")
plt.scatter(x_3, y_3, c='red', s=ps, marker='o', label='3', linewidths=lw, edgecolors="w")
plt.scatter(x_4, y_4, c='blue', s=ps, marker='o', label='4', linewidths=lw, edgecolors="w")
plt.scatter(x_5, y_5, c='black', s=ps, marker='o', label='5', linewidths=lw, edgecolors="w")

sigma    = 10 ## good value for Gaussian kernel
##sigma    = 0.5

n        = 180 

def Gaussian(x1, y1, x2, y2, z2, sigma):
    v_pixel = np.array([x1, y1])
    wtot    = 0.
    wetot   = 0.
    for k in range(len(x2)):
        xk = x2[k]
        yk = y2[k]
        zk = z2[k]
        v_point = np.array([xk, yk])
        delta   = ( v_pixel - v_point )
        d2      = np.dot(delta, delta)
        
        w       = np.exp(-d2/(2*sigma**2))
        #w       = d2**(-sigma)
        
        #print("centre %i at %f %f, E: %f d: %f w: %e" % \
        #                (k, xk, yk, zk, np.sqrt(d2), w))
        
        wetot  += w * zk
        wtot   += w

    return wetot/wtot

Ei = E0
Ej = E1

#definition of min and max of alpha and chi
min_X = np.min(Ei)
max_X = np.max(Ei)
min_Y   = np.min(Ej)
max_Y   = np.max(Ej)

#pixel side
x_space = np.linspace(min_X, max_X, num = n)
y_space = np.linspace(min_Y,   max_Y,   num = n)
E_matrix = np.zeros( (n,n) )


print("scanning pixels over range: ", x_space[0], "  .. ", x_space[-1])
print("                            ", y_space[0], "  .. ", y_space[-1])


#creation of pixel and energy matrix with sum(w*E)/wtot
for i in range(len(x_space)):    
    x = x_space[i]
    for j in range(len(y_space)):
        y = y_space[j]

        #print("pix: %f, %f " % (x,y))
        E_matrix[i,j] = Gaussian(x, y, Ei, Ej, energynm, sigma)  

lim = [min_X, max_X, min_Y, max_Y]
plt.imshow(E_matrix.transpose(),\
           extent = lim,\
           cmap   = 'gist_ncar',\
           origin = "lower",\
           aspect = 'auto')

clb = plt.colorbar()
clb.set_label('Energy error (kcal/mol)', labelpad=-40, y=1.05, rotation=0)
plt.xlabel('e0')
plt.ylabel('e1')
plt.title( 'Energy error in e0-e1 space') 
plt.legend()
plt.savefig("heatmap_energy_error_e0_e1_exact.pdf")
