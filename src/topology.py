import numpy as np

class Topology:

    import numpy as np
import time
import json

class Topology:
    
    def __init__(self, length, width, height, nelx, nely, nelz):
        
        self.length = length
        self.width = width
        self.height = height
        self.nelx = nelx
        self.nely = nely
        self.nelz = nelz
   
        self.readFromCache = False

        self.connections = None
        self.coordinates = None
        self.conelements = None

        self.readCachedContents()
        self.writeCachedContents()
        self.design = np.ones((self.getNumElements(), 1))
        
    def readCachedContents(self):
        
        with open('topologymemory.json', 'r') as topologycache:
            
            topologyinfo = json.load(topologycache)
            key = str(self.nelx) + "x" + str(self.nely) + "x" + str(self.nelz)
            
            if key in topologyinfo:
                self.connections = np.array(topologyinfo[key]["connections"])
                self.coordinates = np.array(topologyinfo[key]["coordinates"])
                self.conelements = np.array(topologyinfo[key]["conelements"])
                print('Topology information read from memory')
                self.readFromCache = True

        topologycache.close()

        print('Leaving read saved memory...')

    def writeCachedContents(self):
        
        if self.readFromCache == False:
            print("Writing topology information to memory...")
            connections = self.getConnections().tolist()
            coordinates = self.getCoordinates().tolist()
            conelements = self.getConnectedElements().tolist()

            dimensions = str(self.nelx) + "x" + str(self.nely) + "x" + str(self.nelz)

            cachedData = {}
            cachedData[dimensions] = {}

            cachedData[dimensions]["connections"] = connections
            cachedData[dimensions]["coordinates"] = coordinates
            cachedData[dimensions]["conelements"] = conelements

            with open('topologymemory.json', 'w+') as topologycache:
                json.dump(cachedData, topologycache)
            
            topologycache.close()

    def getNumElements(self):
        return self.nelx * self.nely * self.nelz
        
    def getElementDimensions(self):
        return (self.nelx, self.nely, self.nelz)
        
    def getNumNodes(self):
        return (self.nelx+1) * (self.nely+1) * (self.nelz+1)

    def getConnections(self):
        
        if self.readFromCache == False:
            connectivity = np.zeros((self.getNumElements(), 8))
        
            for point in range(1, self.getNumElements()+1):
                elx = np.fix((point-1)/(self.nely*self.nelz)) + 1
                ely = np.fix(((point-(elx-1)*self.nely*self.nelz)-1)/self.nelz) + 1
                elz = np.mod(point, self.nelz)
            
                if elz == 0:
                    elz = self.nelz
            
                connectivity[point-1][0] = (self.nely+1)*(self.nelz+1)*(elx-1)+(ely-1)*(self.nelz+1)+elz
                connectivity[point-1][1] = connectivity[point-1][0] + (self.nely+1)*(self.nelz+1)
                connectivity[point-1][2] = connectivity[point-1][1] + (self.nelz+1)
                connectivity[point-1][3] = connectivity[point-1][0] + (self.nelz+1)
                connectivity[point-1][4:8] = (connectivity[point-1][0:4]+1)
            return (connectivity-1).astype(int) 
        else:
            connectivity = self.connections
            
        return connectivity
        
    def getCoordinates(self):
        
        if self.readFromCache == False:

            coordinates = np.zeros((self.getNumConnections(), 3))
        
            dx = self.length / self.nelx
            dy = self.width / self.nely 
            dz = self.height / self.nelz
        
            for point in range(1, self.getNumConnections()+1):
                nx = np.fix((point-1)/((self.nely+1)*(self.nelz+1))) + 1
                ny = np.fix(((point-(nx-1)*(self.nely+1)*(self.nelz+1))-1)/(self.nelz+1)) + 1
                nz = np.mod(point, (self.nelz+1))
            
                if nz == 0:
                    nz = self.nelz + 1
                
                coordinates[point-1] = [(nx-1)*dx, (ny-1)*dy, (nz-1)*dz]
            
            coordinates = coordinates.astype(int)

        else:
            coordinates = self.coordinates

        return coordinates
        