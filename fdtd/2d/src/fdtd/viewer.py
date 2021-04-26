import matplotlib.pyplot as plt
import numpy as np

class View:
    """ Ahora mismo no hace nada"""

    def __init__(self,datos):
        self.data = datos[0]
        self.Nx = len(self.data['values'][0]) # Cantidad de datos en el eje X
        self.X_axis = [self.data['mesh']['origin'][0] + i*self.data['mesh']['steps'][0] \
            for i in range(0,self.Nx)]

        self.Ny = len(self.data['values'][0][0]) # Cantidad de datos en el eje Y
        self.Y_axis = [self.data['mesh']['origin'][1] + i*self.data['mesh']['steps'][1] \
            for i in range(0,self.Ny)]

    def plot(self,time):
        X,Y = np.meshgrid(self.X_axis, self.Y_axis)
        Z = self.data['values'][time]
        
        plt.contour(X, Y, Z)                           
        plt.show()