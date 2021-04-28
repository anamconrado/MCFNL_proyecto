import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import animatplot as amp
import copy

class View:
    """ Allows the visualization of simulation results. """

    def __init__(self, datos):
        self.data = copy.deepcopy(datos[0])
        self.Ntimes = len(self.data['time'])
        self.Nhzx = len(self.data['values'][0]) # Cantidad de datos en el eje X
        self.Nhzy = len(self.data['values'][0][0]) # Cantidad de datos en el eje Y
     

        # Z-component magnetic field (Hz) grid
          # Fields origen
        self.data['mesh']['originh'] = [self.data['mesh']['origin'][0] + self.data['mesh']['steps'][0]/2, \
            self.data['mesh']['origin'][1] + self.data['mesh']['steps'][0]/2]

        self.x_axis = [self.data['mesh']['originh'][0] + i*self.data['mesh']['steps'][0] \
            for i in range(0,self.Nhzx)]
        
        self.y_axis = [self.data['mesh']['originh'][1] + i*self.data['mesh']['steps'][1] \
            for i in range(0,self.Nhzy)]

        """
        # Grid de campo eléctrico Ex:
        self.data['mesh']["origine"] = [self.data['mesh']["origine"][0]\
             + self.data['mesh']['steps']/2, self.data['mesh']["origine"][1]] 
        self.exX_axis = [self.data['mesh']['origin'][0] + i*self.data['mesh']['steps'][0] \
            for i in range(0,self.Nhx)]
        """


    def plot(self, time, fields = "magnetic"):
        """ Plots a snapshot of the input field at the input time:
        Inputs:
        | - time: whatever value between 0 and finalTime. 
        | - fields: must be 'magnetic': plots Hz;'electric': plots de module of the elctric
        |field or 'both': plots both at the same time."""
        X,Y = np.meshgrid(self.x_axis,self.y_axis)
        Z = self.data['values'][time]

        plt.contour(X, Y, Z)                           
        plt.show()

    def Ztimes(self,t):
        return self.data['values'][t]

    def generate_video(self, fields = "magnetic"):
        """ Generates a visualization of the dynamics of the input field and writes a mp4
        video as output:
        Inputs:
        | - fields: must be 'magnetic': plots Hz; 'electric': plots de module of the electric
        |field or 'both': plots both at the same time."""

        # Creating arrays of indices for x-y axis and time.
        x_ind = range(0,len(self.x_axis))
        y_ind = range(0,len(self.y_axis))
        t_ind = range(0,len(self.data["time"]))

        
        # Creating a grid with time and x-y axis data:
        X, Y, _ = np.meshgrid(self.x_axis,self.y_axis,self.data["time"])

        # Creating a figure for visualization
        plt.figure()

        plt.xlabel('x')
        plt.ylabel('y')
        

        if fields == 'magnetic':
            # Se crea una matriz pcolormesh_data con los datos correspondientes a cada punto del grid:
            pcolormesh_data = np.array([np.array([np.array([self.data['values'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])

            plt.title(r'${H_z}$')

            # Animation        
            # now we make our blocks
            pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data,
                                          t_axis=2,vmin = -0.05, vmax = 0.05, cmap='RdBu')
            plt.colorbar(pcolormesh_block.quad)
            timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=30, units='ns')

            # now to contruct the animation
            anim = amp.Animation([pcolormesh_block], timeline)
            anim.controls()

            anim.save('videos/magnetic_Z_field.mp4')
            plt.show()

        elif fields == 'electric':
            pcolormesh_data = np.array([np.array([np.array([self.data['valuese'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])      

            plt.title(r'${ \left | \vec E \right | }$')  

            #Animation   
            # now we make our blocks
            pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data,
                                          t_axis=2,vmin = 0, vmax = 0.05, cmap='RdBu')
            plt.colorbar(pcolormesh_block.quad)
            timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=30, units='ns')

            # now to contruct the animation
            anim = amp.Animation([pcolormesh_block], timeline)
            anim.controls()

            anim.save('videos/electric_field_magnitude.mp4')
            plt.show() 





