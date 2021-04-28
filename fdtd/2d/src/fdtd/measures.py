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
        self.port_inc = self._measures['port_inc']
        self.box_inc = self._mesh.elemIdToBox(self.port_inc["elemId"])
        self.ids_inc = mesh.toIdx(self.box_inc)
        self.port_trans = self._measures['port_trans']
        self.box_trans = self._mesh.elemIdToBox(self.port_trans["elemId"])
        self.ids_trans = mesh.toIdx(self.box_trans)
        self.port_refl = self._measures['port_refl']
        self.box_refl = self._mesh.elemIdToBox(self.port_refl["elemId"])
        self.ids_refl = mesh.toIdx(self.box_refl)

        mu = 1.0
        epsilon = 1.0
        self.eta = np.sqrt(mu/epsilon)
       
        self.Ntimes = len(self._data['time'])
        self.Nhzx = len(self._data['values'][0]) # Cantidad de datos en el eje X
        self.Nhzy = len(self._data['values'][0][0]) # Cantidad de datos en el eje Y

    def I_inc_f(self):
        self.Hz_inc = np.array([np.array([np.array([k for k in j[self.ids_inc[L][Y]:self.ids_inc[U][Y]]]) for j in i[self.ids_inc[L][X]: self.ids_inc[U][X]]]) for i in self._data["values"]])
        self.I_inc = self.eta * np.sum(np.sum(np.power(self.Hz_inc, 2), axis = 1), axis = 1)               
        return self.I_inc

    def I_trans_f(self):
        self.Hz_trans = np.array([np.array([np.array([k for k in j[self.ids_trans[L][Y]:self.ids_trans[U][Y]]]) for j in i[self.ids_trans[L][X]: self.ids_trans[U][X]]]) for i in self._data["values"]])
        self.I_trans = self.eta * np.sum(np.sum(np.power(self.Hz_trans, 2), axis = 1), axis = 1)               
        return self.I_trans

    def I_refl_f(self):
        self.Hz_refl = np.array([np.array([np.array([k for k in j[self.ids_refl[L][Y]:self.ids_refl[U][Y]]]) for j in i[self.ids_refl[L][X]: self.ids_refl[U][X]]]) for i in self._data["values"]])
        self.I_refl = self.eta * np.sum(np.sum(np.power(self.Hz_refl, 2), axis = 1), axis = 1)               
        return self.I_refl 

    def R_f(self):
        R = self.I_refl_f()[1:] / self.I_inc_f()[1:]
        return R

    def T_f(self):
        T = self.I_trans_f()[1:] / self.I_inc_f()[1:]
        return T