import numpy as np
import copy

from fdtd.common import X, Y, L, U

class Measures:
    def __init__(self, mesh, data, measures, material):
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

        self._material = copy.deepcopy(material)
        epsilon = self._material["epsilon"]
        mu = self._material["mu"]
        self.eta = np.sqrt(mu/epsilon)


    """
    def I_inc_f(self):
        self.Hz_inc = np.array([np.array([np.array([k for k in j[self.ids_inc[L][Y]:self.ids_inc[U][Y]]]) for j in i[self.ids_inc[L][X]: self.ids_inc[U][X]]]) for i in self._data["values"]])
        self.I_inc = self.eta * np.sum(np.sum(np.power(self.Hz_inc, 2), axis = 1), axis = 1)               
        return self.I_inc
    """

    def Ports(self,n):
        if n == 0:
            ids =  self.ids_inc
            x_width = float(ids[U][X]-ids[L][X])
        elif n == 2:
            ids = self.ids_refl
            x_width = float(ids[U][X]-ids[L][X])
        elif n == 1:
            ids = self.ids_trans  
            x_width = float(ids[U][X]-ids[L][X])
        else: raise Exception("0, 1 and 2 ports from left to right")
        self.fields = np.array([{"Hz": np.array([np.array([k for k in j[ids[L][Y]:ids[U][Y]]]) for j in self._data["values"][i][ids[L][X]: ids[U][X]]]) ,\
            "Ex": np.array([np.array([k for k in j[ids[L][Y]:ids[U][Y]]]) for j in self._data["valuese_x"][i][ids[L][X]: ids[U][X]]]) ,\
            "Ey": np.array([np.array([k for k in j[ids[L][Y]:ids[U][Y]]]) for j in self._data["valuese_y"][i][ids[L][X]: ids[U][X]]])}\
            for i in range(0,len(self._data["values"]))])
        self.Power = np.array([sum(sum((-1)*i["Ey"]*i["Hz"]))*self._mesh.dy*(1/x_width) for i in self.fields])    
        return self.Power
        
    """
    def I_trans_f(self):
        self.Hz_trans = np.array([np.array([np.array([k for k in j[self.ids_trans[L][Y]:self.ids_trans[U][Y]]]) for j in i[self.ids_trans[L][X]: self.ids_trans[U][X]]]) for i in self._data["values"]])
        self.I_trans = self.eta * np.sum(np.sum(np.power(self.Hz_trans, 2), axis = 1), axis = 1)               
        return self.I_trans
    """

    """
    def I_refl_f(self):
        self.Hz_refl = np.array([np.array([np.array([k for k in j[self.ids_refl[L][Y]:self.ids_refl[U][Y]]]) for j in i[self.ids_refl[L][X]: self.ids_refl[U][X]]]) for i in self._data["values"]])
        self.I_refl = self.eta * np.sum(np.sum(np.power(self.Hz_refl, 2), axis = 1), axis = 1)               
        return self.I_refl 
    """

    def R_f(self):
        # R = self.I_refl_f()[1:] / self.I_inc_f()[1:]
        R = self.Ports(1)[1:]/self.Ports(0)[1:]
        return R

    def T_f(self):
        # T = self.I_trans_f()[1:] / self.I_inc_f()[1:]
        T = self.Ports(2)[1:]/self.Ports(0)[1:]
        return T