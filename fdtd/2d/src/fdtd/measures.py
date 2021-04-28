import numpy as np
import copy

X = 0 # Cartesian indices
Y = 1

L = 0 # Lower
U = 1 # Upper

class Measures:
    def __init__(self, mesh, data, measures):
        self._mesh = copy.deepcopy(mesh)
        self._data = copy.deepcopy(data[0])
        self._measures = copy.deepcopy(measures)
        port_inc = self._measures['port_inc']
        box_inc = self._mesh.elemIdToBox(port_inc["elemId"])
        self.ids_inc = mesh.toIdx(box_inc)
        port_trans = self._measures['port_inc']
        box_trans = self._mesh.elemIdToBox(port_inc["elemId"])
        self.ids_trans = mesh.toIdx(box_inc)
        port_refl = self._measures['port_inc']
        box_refl = self._mesh.elemIdToBox(port_inc["elemId"])
        self.ids_refl = mesh.toIdx(box_inc)
       
        self.Ntimes = len(self._data['time'])
        self.Nhzx = len(self._data['values'][0]) # Cantidad de datos en el eje X
        self.Nhzy = len(self._data['values'][0][0]) # Cantidad de datos en el eje Y

    def I_inc(self):
        mu = 1.0
        epsilon = 1.0
        eta = np.sqrt(mu/epsilon)
        self.I_inc = [self._data['time']]
        self.I_inc.append(np.sum(np.power(np.array(self._data['values'][:] \
                            [self.ids_inc[L][X]:self.ids_inc[U][X]] \
                            [self.ids_inc[L][Y]:self.ids_inc[U][Y]]), 2), axis = 1))                
        return self.I_inc

    def R(self):
        R = self.I_inc()
        return R
