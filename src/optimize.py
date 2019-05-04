import numpy as np
from scipy.sparse import lil_matrix

class Optimizer:

    """ Summary: This class is responsible for optimizing a partiuclar 
        topology for a mean compliance criterion. This class uses
        the bidirectional evolutionary structure optimization (BESO)
        method, but can be extended to work for other sorts of methods """

    def __init__(self, topology, volfrac, force, forceLoc, youngMod, nu=0.4, rmin=1.5, penalization=3):

        self.topology = topology
        self.volfrac = volfrac
        self.force = force
        self.forceLoc = forceLoc
        self.youngMod = youngMod
        self.nu = nu
        self.rmin = rmin
        self.penalization = penalization

        self.volfracit = 0
        self. dc = 0
        self.ndc = 0

        self.MAX_ITERATIONS = 100

    def BESO3D(self):

        #-----------------------------------------------------------------------------
        #                             Optimization Setup
        #-----------------------------------------------------------------------------

        # Get number elements and nodes in system
        numElems = self.topology.getNumElems()
        numNodes = self.topology.getNumNodes()
        connections = self.topology.getConnections()
        coordinates = self.topology.getCoordinates()
        freedofs = self.topology.getFreeDOFS()

        # Displacement and load matrices
        U = lil_matrix((3*numNodes, 1))
        F = np.zeros((3*numNodes, 1))
        design = np.zeros((numElems, 1))

        # Set loads on particular nodes
        F[self.forceLoc-1, 0] = self.force

        # Get the stiffness matrix of entire topology
        Ke, Ve = self.feasolver.getStiffnessMatrix()


        #-----------------------------------------------------------------------------
        #                           Begin Optimization Process
        #-----------------------------------------------------------------------------

        for iteration in range(self.MAX_ITERATIONS):
            
            # Keep track of mean compliance
            oldc = self.dc if iteration > 1 else 0

            # Update the current volume fraction
            self.volfracit = np.maximum(self.volfracit * (1-self.er), self.volfrac)

            # Update the stiffness matrix
            K = self.feasolver.updateK(numElems, numNodes, connections, Ke, design)

            # Solve for nodal displacements, K = UF
            U[freedofs, 0] = np.linalg.solve(K[np.ix_(freedofs, freedofs)], F[freedofs, 0][:, np.newaxis])

            # Get sensitivity values fore ach element in topology
            changeInMC, MC  = self.getSensitivities(numElems, connections, U, design, Ke, Ve)

            # Neighboring node filter
            ncon, ndc, dc = self.neighborFiltering(iteration, numNodes, connections, dc, oldc)

            # Add or substract material
            design = self.addOrDelete(dc, design)

            # Create visualization of structure
            self.visualizer.createVTK(design, iternum=iteration, filetype="vtk", dir=self.filePath)
    



    def getSensitivities(self, numelems, con, U, design, Ke, elemVolume):

        # Elemental sensitivites and change in mean compliance
        ese = np.zeros((numelems, 1))
        dc = np.zeros((numelems, 1))

        # Get all element degrees of freedoms (DOFS) 
        for element in range(numelems):

            elementDOFS = np.array([
                3*con[element][0]-3, 3*con[element][0]-2, 3*con[element][0]-1, 
                3*con[element][1]-3, 3*con[element][1]-2, 3*con[element][1]-1,
                3*con[element][2]-3, 3*con[element][2]-2, 3*con[element][2]-1,
                3*con[element][3]-3, 3*con[element][3]-2, 3*con[element][3]-1,
                3*con[element][4]-3, 3*con[element][4]-2, 3*con[element][4]-1,
                3*con[element][5]-3, 3*con[element][5]-2, 3*con[element][5]-1,
                3*con[element][6]-3, 3*con[element][6]-2, 3*con[element][6]-1,
                3*con[element][7]-3, 3*con[element][7]-2, 3*con[element][7]-1
            ]).astype(int)

            # Calculate mean compliance and sensitivity values
            Ue = U[elementDOFS, 0]
            UeTKe = np.matmul(Ue.T.toarray(), Ke[:,:,element])
            ese[element, 0] = 0.5 * np.matmul(UeTKe, Ue.toarray())[0][0] * design[element]**self.penalization
            dc[element, 0] = ese[element, 0] / elemVolume


        return dc, np.sum(ese)




    def neighborFiltering(self, iteration, numNodes, connection, dc, olddc):

        ncon = np.zeros((numNodes, 8))
        ndc = np.zeros((numNodes, 1))

        # Get node connections
        for node in range(1, numNodes+1):
            row, col = np.where(connection == node)
            row = row[::-1]
            ncon[node-1][0:len(row)] = row+1

        # Get node mean compliance
        for node in range(numNodes):
            a = np.where(ncon[node, :] > 0)
            b = dc[ncon[node, a][0].astype(int)-1, 0]
            ndc[node, 0] = sum(b) / np.maximum(8, len(a))

        oldndc = ndc

        # Calculate average complaince based on filter radius
        for i in range(ndc.shape[0]):
            wi = rmin - b[i]
            alphai = oldndc[a[i][0].astype(int), 0].conj().transpose()
            ndc[i, 0] = sum(wi * alphai) / sum(wi)

        # Assign average values to areas within filter radius
        for i in range(dc.shape[0]):
            dc[i, 0] = np.mean(ndc[connection[i,:].conj().transpose().astype(int)-1, 0])

        # Update change of mean compliance to help convergence (Xie)
        if iteration > 1:
            dc = (dc + olddc)

        return ncon, ndc, dc
        

    def addOrDelete(self, dc, design, Ve, volfrac, domainVol):
        
        limit1, limit2 = min(dc), max(dc)

        # Add or delete elements in design based on threshold
        while((limit2-limit1)/limit2 > 1.0e-5):

            threshold = (limit1+limit2) / 2.0
            design = np.maximum(0.001, np.sign(dc-threshold))

            if np.sum(design * Ve[0]) - volfrac * (domainVol):
                limit1 = threshold
            else:
                limit2 = threshold
        
        return design