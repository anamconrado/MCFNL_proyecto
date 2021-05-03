import numpy as np
import copy

from fdtd.common import X, Y, L, U

class Freq_analysis:
    def __init__(self, measures, data, data_input):
        # Power values
        self._measures = copy.deepcopy(measures)
        
        # Times
        self._times = copy.deepcopy(data[0]["time"])
        
        # Ports positions
        sourceposition = data_input["coordinates"][data_input["elements"][data_input["sources"][0]["elemId"]][0]][0]
        self._ports = {"0": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_inc"]["elemId"]][0]][0])],\
            "1": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_refl"]["elemId"]][0]][0])],\
            "2": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_trans"]["elemId"]][0]][0])],}
        
        # Port times
        self._ports["0"].append([0,10])     # Firsts times at port 0 in ns
        self._ports["1"].append([14,24])   # Second times at port 1 in ns
        self._ports["2"].append([40,50])   # Firsts times at port 3 in ns

        for i in self._ports:
            self._ports[i].append([(self._times[j], self._measures.Ports(int(i))[j])\
                for j in range(0,len(self._times)) if (self._ports[i][1][0] < self._times[j]*10**9 and self._ports[i][1][1] > self._times[j]*10**9)])

        print("hola")

        def Fourier_trans():
            # Calculo de frecuencias
            timestep = self._times[1] - self._times[0]
            for i in self._ports.values():
                i.append((np.fft.fftfreq(len(i[2][0]))/timestep), np.fft.fft(i[2][1]) )
        
