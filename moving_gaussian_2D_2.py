""" Función Gaussiana 2D moviendose en la diagonal usando contour pltos. """
import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import animatplot as amp

# Parameters
twidth = 20
width = 20
high = 20
Nt = 50
Nx = 60
Ny = 60
colormap = colors.Colormap("whiteminima")
colormap.set_under(color = 'w')

def Gaussian(t,x,y,A = 1, sig = 1, c = 1):
    '''
        Computes Gaussian function
    
    '''
    
    x2 = (x - (np.sqrt(2)/2)*c*t)
    y2 = (y - (np.sqrt(2)/2)*c*t)
    return A*np.exp(-(x2**2 + y2**2)/(sig*np.sqrt(2))**2)

t = np.linspace(-twidth/2,twidth/2, Nt)
x = np.linspace(-width/2, width/2, Nx+1)
y = np.linspace(-high/2, high/2, Ny+1)

def Ztimes(t_value):
    return [[Gaussian(t_value,i,j)  for j in y] for i in x]
Z = np.array([Ztimes(i) for i in t])
X, Y, T = np.meshgrid(x, y, t)


pcolormesh_data = Gaussian(T,X,Y)
line_data       = pcolormesh_data[int(Nx/2),:,:] 


# standard matplotlib stuff
# create the different plotting axes
plt.figure()

plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$z=Ae^{b(\vec r-\vec v t)^2}$')

x_ind = range(0,Nx+1)
y_ind = range(0,Ny+1)
t_ind = range(0,Nt)

X_ind, Y_ind, T_ind = np.meshgrid(x_ind, y_ind, t_ind)


pcolormesh_data2 = np.array([np.array([np.array([Z[t][i][j] for t in t_ind]) for i in x_ind]) for j in y_ind])

# animatplot stuff
# now we make our blocks
pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data2,
                                          t_axis=2, vmin=-1, vmax=1)
plt.colorbar(pcolormesh_block.quad)
timeline = amp.Timeline(t, fps=10)

# now to contruct the animation
anim = amp.Animation([pcolormesh_block], timeline)
anim.controls()

anim.save_gif('multiblock')
plt.show()
