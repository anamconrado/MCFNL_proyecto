import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
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
        
        if fields == "magnetic":
            X,Y = np.meshgrid(self.hzx,self.hzy)
            fig = plt.figure()
            ax1 = plt.axes(xlim=(self.hzx[0], self.hzx[-1]), ylim=(self.hzy[0], self.hzy[-1]))
            ax1.set_title("Magnetic field Hz")        
            
            def animate(i):
                Z = self.Ztimes(i)
                ax1.collections = [] 
                cont = plt.contourf(X, Y, Z, levels = 8)
                return cont
            
            anim = animation.FuncAnimation(fig, animate, frames = self.Ntimes, interval = 50)

            anim.save('Magnetic_field.mp4')
