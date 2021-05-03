import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from fdtd.common import X, Y, L, U

def means2(vec):
    lenvec = len(vec)
    if lenvec % 2 == 0:
        vecpar = np.array([vec[i] for i in range(0,lenvec) if i%2 == 0])
        vecimp = np.array([vec[i] for i in range(0,lenvec) if i%2 == 1])
    else: 
        vecpar = np.array([vec[i] for i in range(0,lenvec-1) if i%2 == 0])
        vecimp = np.array([vec[i] for i in range(0,lenvec) if i%2 == 1])       
    return (vecpar + vecimp)/2

def Fourier_trans(measures,data,data_input):
    # Power values
    measures = copy.deepcopy(measures)
    
    # Times
    times = np.array(list(map(lambda x: x*10**9, copy.deepcopy(data[0]["time"]))))
    Ntimes = len(times)

    # Ports positions
    sourceposition = data_input["coordinates"][data_input["elements"][data_input["sources"][0]["elemId"]][0]][0]
    ports = {"0": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_inc"]["elemId"]][0]][0])],\
        "1": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_refl"]["elemId"]][0]][0])],\
        "2": [abs(sourceposition - data_input["coordinates"][data_input["elements"][data_input["measures"]["port_trans"]["elemId"]][0]][0])],}
    
    # Port times
    ports["0"].append([0,10])     # Firsts times at port 0 in ns
    ports["1"].append([15,25])   # Second times at port 1 in ns
    ports["2"].append([41,51])   # Firsts times at port 3 in ns

    for i in ports:
        ports[i].append([[],[]]) 


    for i in ports:
        for j in range(0,Ntimes): 
            if (ports[i][1][0] < times[j] and ports[i][1][1] > times[j]):
                ports[i][2][1].append(abs(measures.Ports(int(i))[j]))
                ports[i][2][0].append(times[j])
            else:
                ports[i][2][1].append(0)
                ports[i][2][0].append(times[j])

    # Nmeasures = len(ports[i][2][0])
    Nmeasures = Ntimes

    # Calculo de frecuencias y transformada de Fourier.
    timestep = times[1] - times[0]
    for i in ports.values():
        frequencies = (fftfreq(Nmeasures)/timestep)[:Nmeasures//2]
        transform = np.abs((fft(i[2][1]))[:Nmeasures//2])
        i.append([frequencies, transform] )

    # Calculo de coeficientes.
    R = means2(ports["2"][3][1])/means2(ports["0"][3][1])
    T = means2(ports["1"][3][1])/means2(ports["0"][3][1])

    return (ports,[R,T])

