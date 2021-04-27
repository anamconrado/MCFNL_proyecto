import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import animatplot as amp
import copy

class View:
    """ Clase que permite visualizar los resultados de la simulación. """

    def __init__(self, datos):
        self.data = copy.deepcopy(datos[0])
        self.Ntimes = len(self.data['time'])
        self.Nhzx = len(self.data['values'][0]) # Cantidad de datos en el eje X
        self.Nhzy = len(self.data['values'][0][0]) # Cantidad de datos en el eje Y
     

        # Grid de campo magnético Hz
          # Origen de campo magnético
        self.data['mesh']['originh'] = [self.data['mesh']['origin'][0] + self.data['mesh']['steps'][0]/2, \
            self.data['mesh']['origin'][1] + self.data['mesh']['steps'][0]/2]

        self.hzx = [self.data['mesh']['originh'][0] + i*self.data['mesh']['steps'][0] \
            for i in range(0,self.Nhzx)]
        
        self.hzy = [self.data['mesh']['originh'][1] + i*self.data['mesh']['steps'][1] \
            for i in range(0,self.Nhzy)]

        """
        # Grid de campo eléctrico Ex:
        self.data['mesh']["origine"] = [self.data['mesh']["origine"][0]\
             + self.data['mesh']['steps']/2, self.data['mesh']["origine"][1]] 
        self.exX_axis = [self.data['mesh']['origin'][0] + i*self.data['mesh']['steps'][0] \
            for i in range(0,self.Nhx)]
        """


    def plot(self, time, fields = "magnetic"):
        X,Y = np.meshgrid(self.hzx,self.hzy)
        Z = self.data['values'][time]

        plt.contour(X, Y, Z)                           
        plt.show()

    def Ztimes(self,t):
        return self.data['values'][t]

    def generate_video(self, fields = "magnetic"):
        # Se crean arrays de indices para los ejes x e y y el tiempo: 
        x_ind = range(0,len(self.hzx))
        y_ind = range(0,len(self.hzy))
        t_ind = range(0,len(self.data["time"]))

        # Se crea un grid con los datos del tiempo y los ejes X e Y:
        X, Y, _ = np.meshgrid(self.hzx,self.hzy,self.data["time"])

        # Se crea una matriz pcolormesh_data con los datos correspondientes a cada punto del grid:
        pcolormesh_data = np.array([np.array([np.array([self.data['values'][t][i][j]\
             for t in t_ind]) for i in x_ind]) for j in y_ind])

        # Se animan los datos:

        plt.figure()

        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(r'${H_z}$')


        # now we make our blocks
        pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data,
                                          t_axis=2,vmin = -0.05, vmax = 0.05)
        plt.colorbar(pcolormesh_block.quad)
        timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=10, units='ns')

        # now to contruct the animation
        anim = amp.Animation([pcolormesh_block], timeline)
        anim.controls()

        anim.save_gif('magnetic_field')
        plt.show()