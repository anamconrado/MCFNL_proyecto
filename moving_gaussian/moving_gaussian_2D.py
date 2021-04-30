""" Función Gaussiana 2D moviendose en la diagonal usando contour pltos. """

import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Parameters
twidth = 20
width = 20
high = 20
Nt = 50
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
x = np.linspace(-width/2, width/2, 60)
y = np.linspace(-high/2, high/2, 60)

X, Y = np.meshgrid(x,y)

def Ztimes(t_value):
    return [[Gaussian(t_value,i,j)  for j in x] for i in y]

Z = Ztimes(t)

x0 = np.array([])
y0 = np.array([])
X0, Y0 = np.meshgrid(x0,y0)
Z0 = np.zeros((len(x0),len(y0)))

fig = plt.figure()
ax1 = plt.axes(xlim=(x[0], x[-1]), ylim=(y[0], y[-1]))
ax1.set_title("Moving Gaussian")

ax1.grid()

def animate(i):
    Z = Ztimes(t[i])
    ax1.collections = [] 
    cont = plt.contourf(X, Y, Z, levels = 5, colors = ["#FFFFFF","#CCCC00","#FF9933","#FF8000","#FF3333"])

    return cont

anim = animation.FuncAnimation(fig, animate, frames=Nt, interval = 100)

anim.save('animation.mp4')

