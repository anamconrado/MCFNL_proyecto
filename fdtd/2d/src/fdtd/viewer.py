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


    def plots(self, measures):
        """ Plots the input of a port defined in the measures module:
        Inputs:
        | - measures: Object of Measures, give access to class functions as Ports.
        | - port: Must be 0,1 or 2: indicates the port whose data will be plotted.
        Output:
        | - Plot of input port data at all times.
        """
        fig, (ax1, ax2, ax3) = plt.subplots(3,1)
        ports = {"0":[ax1,0] ,"1":[ax2,1] ,"2":[ax3,2] }
        for i in ports:
            ports[i][0].plot(list(map(lambda i: i*(10**9), self.data['time'])), measures.Ports(ports[i][1]))
            ports[i][0].set_title(f"Port {ports[i][1]}")
            ports[i][0].set_ylabel("Power per length")
        ax3.set_xlabel("Time (ns")
        fig.subplots_adjust(left=None, bottom=0.1, right=None, top=0.85, wspace=None, hspace=0.5)
        fig.savefig("Puertos.png")

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
        fig, (ax1, ax2, ax3) = plt.subplots(3,1)

        fields = {"Ex": ['valuese_x',ax1,[],[]], "Ey": ['valuese_y',ax2,[],[]], "Hz":['values',ax3,[],[]]}
        for i in fields.values():
            i[1].set_ylabel('y')
            
        fields["Hz"][1].set_xlabel('x')
        fields["Ex"][1].set_title(r'$ {E_x} $')
        fields["Ey"][1].set_title(r'$ {E_y} $')
        fields["Hz"][1].set_title (r'$ H_z $')


        for i in fields:
            field = fields[i][0]
            fields[i][2] = np.array([np.array([np.array([self.data[field][t][i][j]\
                for t in t_ind]) for i in x_ind]) for j in y_ind]) 

        dicmaxmins = {"Ex": ['valuese_x', [], []],\
                "Ey": ['valuese_y', [], []],\
                "Hz": ['values', [], []]}
        
        t0 = 0
        tf = 100

        for i in dicmaxmins.values():
            for time in self.data[i[0]][t0:tf]: 
                (i[1]).append(max([max(j) for j in time]))
                (i[2]).append(min([min(j) for j in time]))

        for i in dicmaxmins:
            dicmaxmins[i][1] = max(dicmaxmins[i][1])
            dicmaxmins[i][2] = min(dicmaxmins[i][2])

        fig.suptitle(r'${   {E_x}  \ & \  {E_y} \ & \ H_z }$')  
        fig.subplots_adjust(left=None, bottom=0.1, right=None, top=0.85, wspace=None, hspace=0.5)
        #Animation   
        # now we make our blocks
        for i in fields:
            fields[i][3] = amp.blocks.Pcolormesh(X[:,:,0], Y[:,:,0], fields[i][2],
                ax=fields[i][1], t_axis=2,vmin = (dicmaxmins[i][2])*0.6, vmax =  (dicmaxmins[i][1])*0.6)                              

        for i in fields:
            plt.colorbar(fields[i][3].quad, ax = fields[i][1])

        timeline = amp.Timeline([i*(10**9) for i in self.data["time"]], fps=10, units='ns')

        # now to contruct the animation
        anim = amp.Animation([fields["Ex"][3],fields["Ey"][3],fields["Hz"][3],], timeline)
        anim.controls()

        # Change if windows.
        #anim.save_gif('2d/videos/allfselds')
        anim.save('2d/videos/allfields.avi')
        plt.show()                         