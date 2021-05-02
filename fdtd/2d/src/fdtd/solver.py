import math
import numpy as np
import scipy.constants as sp
import copy
import time

from fdtd.common import X, Y, L, U

def gaussian(x, delay, spread):
    return np.exp( - ((x-delay)**2 / (2*spread**2)) )

def step(x, xlim):
    return x<xlim

def subsId(id):
    if id is None:
        return -1
    else:
        return id-1

class Solver:
    
    class Fields: 
        def __init__(self, ex, ey, hz):
            self.ex = ex
            self.ey = ey
            self.hz = hz
        
        def get(self):
            return (self.ex, self.ey, self.hz)

    __timeStepPrint = 5000

    def __init__(self, mesh, options, probes, sources, material):
        self.options = options
        
        self._mesh = copy.deepcopy(mesh)
        (dX, dY) = self._mesh.steps()
        self.posEx = [np.delete(mesh.pos[X] - dX/2, 0), mesh.pos[Y]]
        self.posEy = [mesh.pos[X], np.delete(mesh.pos[Y] - dX/2, 0)]
        self.posHz = [np.delete(mesh.pos[X] - dX/2, 0), np.delete(mesh.pos[Y] - dX/2, 0)] 

        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIdx(box)
            Nxy = abs(ids[Y] - ids[X])
            p["mesh"] = {"origin": box[L], "steps": (dX, dY), "posEx": self.posEx, \
                            "posEy": self.posEy, "posHz": self.posHz}
            p["indices"] = ids
            p["time"]   = [0.0]
            p["values"] = [np.zeros((Nxy[X], Nxy[Y]))]
            p["valuese_mod"] = [np.zeros((Nxy[X], Nxy[Y]))]
            p["valuese_x"] = [np.zeros((Nxy[X], Nxy[Y]))]
            p["valuese_y"] = [np.zeros((Nxy[X], Nxy[Y]))]

        # for initial in self._initialCond:
        #     if initial["type"] == "gaussian":
        #         position=self._mesh.pos
        #         values=Solver.movingGaussian(position, 0, \
        #             sp.speed_of_light,initial["peakPosition"],\
        #             initial["gaussianAmplitude"], \
        #             initial["gaussianSpread"] )  
        #         p["values"]= [values[ids[0]:ids[1]]]
        #     else:
        #         raise ValueError(\
        #         "Invalid initial condition type: " + initial["type"] )

        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIdx(box)
            source["index"] = ids

        self._material = copy.deepcopy(material)
        
        self.old = self.Fields(
            ex = np.zeros( (mesh.pos[X].size-1, mesh.pos[Y].size  ) ),
            ey = np.zeros( (mesh.pos[X].size,   mesh.pos[Y].size-1) ),
            hz = np.zeros( (mesh.pos[X].size-1, mesh.pos[Y].size-1) ) )


    def _dt(self):
        return self.options["cfl"] * min(self._mesh.steps()) / math.sqrt(2.0)  

    def timeStep(self):
        return self._dt() / sp.speed_of_light

    def getProbes(self):
        res = self._probes
        return res

    # ======================= UPDATE E =============================
    def _updateE(self, t, dt, overFields = None):
        eNew = (np.zeros( self.old.ex.shape ),
                np.zeros( self.old.ey.shape ) )
        (ex, ey, h) = self.old.get()
        e = (ex, ey)

        (dX, dY) = self._mesh.steps()
        A = dX * dY
        eNew[X][:,1:-1] = e[X][:,1:-1] + dt/A*dX * (h[:,1:] - h[:,:-1])
        eNew[Y][1:-1,:] = e[Y][1:-1,:] - dt/A*dY * (h[1:,:] - h[:-1,:])

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian": # No hace nada
                    break

                elif magnitude["type"] == "TMgauss":
                    mode = source["mode"]
                    freq = source["frequency"]

                    c0 = sp.speed_of_light
                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    id = source["index"]
                    intens = magnitude["sinIntensity"]

                    epsilon = self._material["epsilon"]
                    mu = self._material["mu"]

                    lon_y = (id[U][Y] - id[L][Y]) * dY
                    kc = mode*np.pi/lon_y
                    beta = np.sqrt(freq**2*mu*epsilon - kc**2)

                    (xEx, yEx) = self._mesh.IdxToPos(id,self.posEx)
                    xEx = xEx[...,None]
                    yEx = yEx[None,...]
                    tEx = t - dt/2
                    eNew[X][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     np.matmul((np.cos(freq*tEx)* np.cos(beta*xEx) - np.sin(freq*tEx) * np.sin(beta*xEx)), \
                     np.sin(kc * yEx)) * intens * gaussian(tEx, delay, spread)

                    (xEy, yEy) = self._mesh.IdxToPos(id,self.posEy)
                    xEy = xEy[...,None]
                    yEy = yEy[None,...]
                    tEy = t - dt/2
                    eNew[Y][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     np.matmul((np.sin(freq*tEy) * np.cos(beta*xEy) + np.sin(beta*xEy) * np.cos(freq*tEy)),
                     np.cos(kc *yEy)) * intens * beta/kc * gaussian(tEy, delay, spread) 

                elif magnitude["type"] == "TMstep":
                    mode = source["mode"]
                    freq = source["frequency"]

                    id = source["index"]
                    intens = magnitude["sinIntensity"]
                    tlim = magnitude["stepTimeLimit"]

                    epsilon = self._material["epsilon"]
                    mu = self._material["mu"]

                    lon_y = (id[U][Y] - id[L][Y]) * dY
                    kc = mode*np.pi/lon_y
                    beta = np.sqrt(freq**2*mu*epsilon - kc**2)

                    (xEx, yEx) = self._mesh.IdxToPos(id,self.posEx)
                    xEx = xEx[...,None]
                    yEx = yEx[None,...]
                    tEx = t + dt/2
                    eNew[X][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     np.matmul((np.cos(freq*tEx)* np.cos(beta*xEx) - np.sin(freq*tEx) * np.sin(beta*xEx)), \
                     np.sin(kc * yEx)) *intens * step(tEx, tlim * dt)

                    (xEy, yEy) = self._mesh.IdxToPos(id,self.posEy)
                    xEy = xEy[...,None]
                    yEy = yEy[None,...]
                    tEy = t + dt/2
                    eNew[Y][id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     np.matmul((np.sin(freq*tEy) * np.cos(beta*xEy) + np.sin(beta*xEy) * np.cos(freq*tEy)),
                     np.cos(kc *yEy)) * intens * beta/kc * step(tEy, tlim * dt)
                     
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
            else:
                raise ValueError("Invalid source type: " + source["type"])

        # Boundary conditions
        for bound in self._mesh.bounds:
            xy = bound.orientation()
            (lx, ux) = (bound.arrayIdx(L,X), \
                        bound.arrayIdx(U,X))
            (ly, uy) = (bound.arrayIdx(L,Y), \
                        bound.arrayIdx(U,Y))

            if isinstance(bound, self._mesh.BoundPEC):
                 eNew[xy][lx:ux,ly:uy] = 0.0
            

#            elif isinstance(bound, self._mesh.BoundMUR):
#                c0 = sp.speed_of_light
#                dx = self._mesh.steps()
#                if ind[0] == 0:
#                    eNew[X][0,:] = e[X][1,:] + \
#                       (c0 * dt - dx)/(c0 * dt + dx)  * (eNew[X][1,:]-e[X][0,:])
#                    eNew[Y][0,:] = e[Y][1,:] + \
#                       (c0 * dt - dx)/(c0 * dt + dx)  * (eNew[Y][1,:]-e[Y][0,:]) 
#                elif ind[0] == -1:
#                    eNew[X][ind,:] = e[X][ind-1,:] + \
#                       (c0 * dt - dx)/(c0 * dt + dx)  * (eNew[X][ind-1,:]-e[X][ind,:])
#                    eNew[Y][ind,:] = e[Y][ind-1,:] + \
#                       (c0 * dt - dx)/(c0 * dt + dx)  * (eNew[Y][ind-1,:]-e[Y][ind,:])
            
            else:
                raise ValueError("Unrecognized boundary type")


        # Subgridding and updating
        e[X][:] = eNew[X][:]
        e[Y][:] = eNew[Y][:]  

    # ======================= UPDATE H =============================
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.hz.shape )
        (ex, ey, h) = self.old.get()
        
        (dX, dY) = self._mesh.steps()
        A = dX * dY
              
        hNew[:,:] = h[:,:] \
                     - dt/A * dY * ey[1:,  :] \
                     + dt/A * dX * ex[ :, 1:] \
                     + dt/A * dY * ey[:-1,   :] \
                     - dt/A * dX * ex[  :, :-1]
        
        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    c0 = sp.speed_of_light
                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    id = source["index"]
                    hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     gaussian(t, delay, spread)*dt

                elif magnitude["type"] == "TMgauss":
                    mode = source["mode"]
                    freq = source["frequency"]

                    c0 = sp.speed_of_light
                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    id = source["index"]
                    intens = magnitude["sinIntensity"]

                    epsilon = self._material["epsilon"]
                    mu = self._material["mu"]

                    lon_y = (id[U][Y] - id[L][Y]) * dY
                    kc = mode*np.pi/lon_y
                    beta = np.sqrt(freq**2*mu*epsilon - kc**2)

                    (xHz, yHz) = self._mesh.IdxToPos(id,self.posHz)
                    xHz = xHz[...,None]
                    yHz = yHz[None,...]
                    hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     np.matmul((-np.sin(freq*t)*np.cos(beta*xHz) - np.sin(beta*xHz)*np.cos(freq*t)), \
                     np.cos(kc * yHz)) * intens * freq*epsilon/kc * gaussian(t, delay, spread)

                elif magnitude["type"] == "TMstep":
                    mode = source["mode"]
                    freq = source["frequency"]

                    id = source["index"]
                    intens = magnitude["sinIntensity"]
                    tlim = magnitude["stepTimeLimit"]

                    epsilon = self._material["epsilon"]
                    mu = self._material["mu"]

                    lon_y = (id[U][Y] - id[L][Y]) * dY
                    kc = mode*np.pi/lon_y
                    beta = np.sqrt(freq**2*mu*epsilon - kc**2)

                    (xHz, yHz) = self._mesh.IdxToPos(id,self.posHz)
                    xHz = xHz[...,None]
                    yHz = yHz[None,...]
                    hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                     np.matmul((-np.sin(freq*t)*np.cos(beta*xHz) - np.sin(beta*xHz)*np.cos(freq*t)), \
                     np.cos(kc * yHz)) * intens * freq*epsilon/kc * step(t, tlim * dt)
                     
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
            else:
                raise ValueError("Invalid source type: " + source["type"])
        
        h[:] = hNew[:]

    def _updateProbes(self, t):
        for p in self._probes:
            dimensionalTime = t/sp.speed_of_light
            writeStep = "samplingPeriod" in p \
                and (dimensionalTime/p["samplingPeriod"] >= len(p["time"]))
            writeStep = writeStep or "samplingPeriod" not in p
            if writeStep:
                p["time"].append(dimensionalTime)
                idx = p["indices"]
                values = np.zeros(tuple(idx[U]-idx[L]))
                values[:,:] = \
                    self.old.hz[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ]
                
                
                # Electric field at magnetic field positions
                
                # Ex values without first raw    
                valuesexup = \
                    np.array(self.old.ex[ idx[L][X]:idx[U][X], (idx[L][Y]+1):(idx[U][Y]+1) ])
                # Ex values without last raw
                valuesexdown = \
                    np.array(self.old.ex[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ])
                

                # Mean values, Ex same position as Hz
                valuesex = (valuesexup + valuesexdown)/2
                
                # Electric field at magnetic field positions
                
                # Ey values without first column
                valueseyright = \
                    self.old.ey[ (idx[L][X]+1):(idx[U][X]+1), idx[L][Y]:idx[U][Y] ]
                # Ey values without last comlumn
                valueseyleft = \
                    self.old.ey[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ]
                

                # Mean values, Ey same position as Hz
                valuesey = (valueseyright + valueseyleft)/2

                p["values"].append(values)
                p["valuese_mod"].append(np.array([list(map(lambda x,y: np.sqrt(x**2 +y**2), valuesex[i], valuesey[i])) for i in range(0,len(valuesex))]))
                p["valuese_x"].append(valuesex)
                p["valuese_y"].append(valuesey)


    def solve(self, dimensionalFinalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        numberOfTimeSteps = \
            int(dimensionalFinalTime * sp.speed_of_light / dt)
        for n in range(numberOfTimeSteps):
        
            self._updateE(t, dt, self.old)
            t += dt/2.0

            self._updateH(t, dt)
            t += dt/2.0

            self._updateProbes(t)
    
            if n % self.__timeStepPrint == 0 or n+1 == numberOfTimeSteps:
                remaining = (time.time() - tic) * \
                    (numberOfTimeSteps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("    Step: %6d of %6d. Remaining: %2.0f:%02.0f"% (n, \
                    numberOfTimeSteps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))

    @staticmethod
    def movingGaussian(x,y,t,c,center,A,spread):
        return A*np.exp(-(((x-center)-c*t)**2 /(2*spread**2)))