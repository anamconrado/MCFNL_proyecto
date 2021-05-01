import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import animatplot as amp
import copy

def Ztimes(self,t):
    return self.data['values'][t]


class View:
    """ Allows the visualization of simulation results. """

    def __init__(self, datos, coeff):
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


    def plots(self, port, measures):
        """ Plots a snapshot of the input field at the input time:
        Inputs:
        | - time: whatever value between 0 and finalTime. 
        | - fields: must be 'magnetic': plots Hz;'electric': plots de module of the elctric
        |field or 'both': plots both at the same time."""
        plt.plot(list(map(lambda i: i*(10**9), self.data['time'])), measures.Ports(port))
        plt.savefig("Puerto {}.png".format(port))

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


        if fields == 'magnetic':
            plt.figure()

            plt.xlabel('x')
            plt.ylabel('y')
        

            # pcolormesh_data matrix with magnited field values.
            pcolormesh_data = np.array([np.array([np.array([self.data['values'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])

            plt.title(r'${H_z}$')

            # Animation        
            # now we make our blocks
            pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data,
                                          t_axis=2,vmin = -0.05, vmax = 0.05)
            plt.colorbar(pcolormesh_block.quad)
            timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=30, units='ns')

            # now to contruct the animation
            anim = amp.Animation([pcolormesh_block], timeline)
            anim.controls()

            anim.save('videos/magnetic_z_field.mp4')
            plt.show()

        elif fields == 'electric':
            plt.figure()
            plt.xlabel('x')
            plt.ylabel('y')

            pcolormesh_data = np.array([np.array([np.array([self.data['valuese_x'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])      

            plt.title(r'${ \vec {E_x} }$')  

            #Animation   
            # now we make our blocks
            pcolormesh_block = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data,
                                          t_axis=2,vmin = -0.05, vmax = 0.05)
            plt.colorbar(pcolormesh_block.quad)
            timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=30, units='ns')

            # now to contruct the animation
            anim = amp.Animation([pcolormesh_block], timeline)
            anim.controls()

            anim.save_gif('videos/electric_field_magnitude')
            plt.show() 
        
        elif fields == 'all':
            fig, (ax1, ax2, ax3) = plt.subplots(3,1)

            for i in [ax1,ax2,ax3]:
                i.set_ylabel('y')
                
            ax3.set_xlabel('x')
            ax3.set_title (r'$ H_z $')
            ax2.set_title(r'$ {E_y} $')
            ax1.set_title(r'$ {E_x} $')

            pcolormesh_data_ex = np.array([np.array([np.array([self.data['valuese_x'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])   
            pcolormesh_data_ey = np.array([np.array([np.array([self.data['valuese_y'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])     
            pcolormesh_data_m = np.array([np.array([np.array([self.data['values'][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind])

            maxe_x = []
            mine_x = []

            maxe_y = []
            mine_y = []
            
            max_m = []
            min_m = []
            t0 = 0
            tf = 100
            for time in self.data['valuese_y'][t0:tf]:
                maxe_y.append(max([max(i) for i in time]))
                mine_y.append(min([min(i) for i in time]))   

            for time in self.data['valuese_x'][t0:tf]:
                maxe_x.append(max([max(i) for i in time]))
                mine_x.append(min([min(i) for i in time]))   

            for time in self.data['values'][t0:tf]:
                max_m.append(max([max(i) for i in time]))
                min_m.append(min([min(i) for i in time]))

            maxmin_e_x = [max(maxe_x),min(mine_x)]
            maxmin_e_y = [max(maxe_y),min(mine_y)]
            maxmin_m = [max(max_m),min(min_m)]

            fig.suptitle(r'${   {E_x}  \ & \  {E_y} \ & \ H_z }$')  
            fig.subplots_adjust(left=None, bottom=0.1, right=None, top=0.85, wspace=None, hspace=0.5)
            #Animation   
            # now we make our blocks
            pcolormesh_block_ex = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data_ex,
                                          ax=ax1, t_axis=2,vmin = maxmin_e_x[1]*0.6, vmax = maxmin_e_x[0]*0.6)
            pcolormesh_block_ey = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data_ey,
                                          ax=ax2, t_axis=2,vmin = maxmin_e_y[1]*0.6, vmax = maxmin_e_y[0]*0.6)
            pcolormesh_block_m = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], pcolormesh_data_m,
                                          ax=ax3, t_axis=2,vmin = maxmin_m[1]*0.6, vmax = maxmin_m[0]*0.6)                                    

            plt.colorbar(pcolormesh_block_ex.quad, ax = ax1)
            plt.colorbar(pcolormesh_block_ey.quad, ax = ax2)
            plt.colorbar(pcolormesh_block_m.quad, ax = ax3)
            timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=10, units='ns')

            # now to contruct the animation
            anim = amp.Animation([pcolormesh_block_ex,pcolormesh_block_ey,pcolormesh_block_m,], timeline)
            anim.controls()

            # Change if windows.
            anim.save_gif('videos/allfields')
            plt.show()                         
        else: raise Exception("Input must be 'magnetic', 'electric' or 'both'")





