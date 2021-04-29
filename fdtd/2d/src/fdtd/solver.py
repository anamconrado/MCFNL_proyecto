import math
import numpy as np
import scipy.constants as sp
import copy
import time

X = 0 # Cartesian indices
Y = 1

L = 0 # Lower
U = 1 # Upper

def gaussian(x, delay, spread):
    return np.exp( - ((x-delay)**2 / (2*spread**2)) )

def sins(x, intens, long):
    return intens*np.sin(np.pi*x/long)

def coss(x, intens, long):
    return intens*np.cos(np.pi*x/long)

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

    def __init__(self, mesh, options, probes, sources):
        self.options = options
        
        self._mesh = copy.deepcopy(mesh)

        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIdx(box)
            Nxy = abs(ids[Y] - ids[X])
            p["mesh"] = {"origin": box[L], "steps": abs(box[U]-box[L]) / Nxy}
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
                    c0 = sp.speed_of_light
                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    id = source["index"]
                    intens = magnitude["sinIntensity"]
                    lon_y = id[U][Y] - id[L][Y]
                    middle_x = int((id[U][X] - id[L][X])/2 + id[L][X])
                    # eNew[X][middle_x, id[L][Y]:id[U][Y]] += \
                    #  sins(np.arange(lon_y), intens, lon_y) * gaussian(t, delay, spread) 
                    # eNew[Y][middle_x, id[L][Y]:id[U][Y]] += \
                    #  coss(np.arange(lon_y), intens, lon_y) * gaussian(t, delay, spread) 

                elif magnitude["type"] == "TMstep":
                    id = source["index"]
                    intens = magnitude["sinIntensity"]
                    tlim = magnitude["stepTimeLimit"]
                    lon_y = id[U][Y] - id[L][Y]
                    middle_x = int((id[U][X] - id[L][X])/2 + id[L][X])
                    #eNew[X][middle_x, id[L][Y]:id[U][Y]] += \
                    # sins(np.arange(lon_y), intens, lon_y) * step(t, tlim * dt) 
                    #eNew[Y][middle_x, id[L][Y]:id[U][Y]] += \
                    # coss(np.arange(lon_y), intens, lon_y) * step(t, tlim * dt) 
                     
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
            #print(lx, ux, ly, ux)

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
                    c0 = sp.speed_of_light
                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    id = source["index"]
                    intens = magnitude["sinIntensity"]
                    lon_y = id[U][Y] - id[L][Y]
                    middle_x = int((id[U][X] - id[L][X])/2 + id[L][X])
                    hNew[middle_x, id[L][Y]:id[U][Y]] += \
                     coss(np.arange(lon_y), intens, lon_y) * gaussian(t, delay, spread) 

                elif magnitude["type"] == "TMstep":
                    id = source["index"]
                    intens = magnitude["sinIntensity"]
                    tlim = magnitude["stepTimeLimit"]
                    lon_y = id[U][Y] - id[L][Y]
                    middle_x = int((id[U][X] - id[L][X])/2 + id[L][X])
                    hNew[middle_x, id[L][Y]:id[U][Y]] += \
                     coss(np.arange(lon_y), intens, lon_y) * step(t, tlim * dt * 10)  
                     
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
                
                """
                # Ex values without first raw    
                valuesexup = \
                    self.old.ex[ idx[L][X]:idx[U][X], (idx[L][Y]+1):(idx[U][Y]+1) ]
                # Ex values without last raw
                valuesexdown = \
                    self.old.ex[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ]
                """

                # Mean values, Ex same position as Hz
                valuesex = np.array([list(map(lambda x, y: (x+y)/2,\
                self.old.ex[ idx[L][X]:idx[U][X], (idx[L][Y]+1):(idx[U][Y]+1) ][i],\
                self.old.ex[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ][i]))\
                for i in range(0,self.old.ex.shape[0])])
                
                """
                # Ey values without first column
                valueseyright = \
                    self.old.ey[ (idx[L][X]+1):(idx[U][X]+1), idx[L][Y]:idx[U][Y] ]
                # Ey values without last comlumn
                valueseyleft = \
                    self.old.ey[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ]
                """

                # Mean values, Ey same position as Hz
                valuesey = np.array([list(map(lambda x, y: (x+y)/2,\
                self.old.ey[ (idx[L][X]+1):(idx[U][X]+1), idx[L][Y]:idx[U][Y] ][i],\
                self.old.ey[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ][i]))\
                for i in range(0,self.old.ey.shape[0]-1)])

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