from elasticity import Elasticity
from pyhull.convex_hull import ConvexHull
from multiprocessing.dummy import Pool as ThreadPool
import scipy.sparse as sp
from threading import Thread
import numpy as np
import time

class FEA:

    def __init__(self, topology, ym, pr):
        
        self.topology = topology
        self.ym, self.pr = ym, pr
        self.con = self.topology.getConnections()
        self.coord = self.topology.getCoordinates()
        self.elasticitymatrix = Elasticity(self.ym, self.pr).getElasticityMatrix()
        self.edofs = []
        self.calcEDOFS()

    def hexElement(self, con, coord, element):

        # Gaussian integration points of hexahedral cube
        g = 1 / np.sqrt(3)

        GP = np.array([
            [-g, -g, -g],
            [g, -g, -g],
            [g, g, -g],
            [-g, g, -g],
            [-g, -g, g],
            [g, -g, g],
            [g, g, g],
            [-g, g, g]])

        i, j, k, l, m, n, o, p = con[element][0:8].astype(int)
        
        # Physical coordinates
        x1, y1, z1 = coord[i-1][0:3].astype(int)
        x2, y2, z2 = coord[j-1][0:3].astype(int)
        x3, y3, z3 = coord[k-1][0:3].astype(int)
        x4, y4, z4 = coord[l-1][0:3].astype(int)
        x5, y5, z5 = coord[m-1][0:3].astype(int)
        x6, y6, z6 = coord[n-1][0:3].astype(int)
        x7, y7, z7 = coord[o-1][0:3].astype(int)
        x8, y8, z8 = coord[p-1][0:3].astype(int)

        ke = np.zeros((24, 24))

        for point in range(8):

            # Integration points
            alpha = round(GP[point][0], 5)
            beta = round(GP[point][1], 5)
            gamma = round(GP[point][2], 5)

            # Lagrangian shape functions for hexahedral cube
            N1 = (1.0/8.0) * (1.0-alpha) * (1.0-beta) * (1.0-gamma)
            N2 = (1.0/8.0) * (1.0+alpha) * (1.0-beta) * (1.0-gamma)
            N3 = (1.0/8.0) * (1.0+alpha) * (1.0+beta) * (1.0-gamma)
            N4 = (1.0/8.0) * (1.0-alpha) * (1.0+beta) * (1.0-gamma)
            N5 = (1.0/8.0) * (1.0-alpha) * (1.0-beta) * (1.0+gamma)
            N6 = (1.0/8.0) * (1.0+alpha) * (1.0-beta) * (1.0+gamma)
            N7 = (1.0/8.0) * (1.0+alpha) * (1.0+beta) * (1.0+gamma)
            N8 = (1.0/8.0) * (1.0-alpha) * (1.0+beta) * (1.0+gamma)

            # Derivatives of Lagrangian shape functions
            dN1dalpha = -(1.0/8.0) * (beta-1.0) * (gamma-1.0)
            dN1dbeta = -(1.0/8.0) * (alpha-1.0) * (gamma-1.0)
            dN1dgamma = -(1.0/8.0) * (alpha-1.0) * (beta-1.0)
            dN2dalpha = (1.0/8.0) * (1.0-beta) * (1.0-gamma)
            dN2dbeta = (1.0/8.0) * (alpha+1.0) * (gamma-1.0)
            dN2dgamma = (1.0/8.0) * (alpha+1.0) * (beta-1.0)
            dN3dalpha = (1.0/8.0) * (1.0+beta) * (1.0-gamma)
            dN3dbeta = (1.0/8.0) * (1.0+alpha) * (1.0-gamma)
            dN3dgamma = -(1.0/8.0) * (alpha+1.0) * (beta+1.0)
            dN4dalpha = (1.0/8.0) * (beta+1.0) * (gamma-1.0)
            dN4dbeta = (1.0/8.0) * (alpha-1.0) * (gamma-1.0)
            dN4dgamma = (1.0/8.0) * (alpha-1.0) * (beta+1.0)
            dN5dalpha = (1.0/8.0) * (beta-1.0) * (gamma+1.0)
            dN5dbeta = (1.0/8.0) * (alpha-1.0) * (gamma+1.0)
            dN5dgamma = (1.0/8.0) * (alpha-1.0) * (beta-1.0)
            dN6dalpha = (1.0/8.0) * (1.0-beta) * (gamma+1.0)
            dN6dbeta = -(1.0/8.0) * (alpha+1.0) * (gamma+1.0)
            dN6dgamma = (1.0/8.0) * (alpha+1.0) * (1.0-beta)
            dN7dalpha = (1.0/8.0) * (beta+1.0) * (gamma+1.0)
            dN7dbeta = (1.0/8.0) * (alpha+1.0) * (gamma+1.0)
            dN7dgamma = (1.0/8.0) * (alpha+1.0) * (beta+1.0)
            dN8dalpha = -(1.0/8.0) * (beta+1.0) * (gamma+1.0)
            dN8dbeta = (1.0/8.0) * (1.0-alpha) * (gamma+1.0)
            dN8dgamma = (1.0/8.0) * (1.0-alpha) * (beta+1.0)

            # Interpolations of shape function
            xalpha = dN1dalpha*x1 + dN2dalpha*x2 + dN3dalpha*x3 + dN4dalpha*x4 + dN5dalpha*x5 + dN6dalpha*x6 + dN7dalpha*x7 + dN8dalpha*x8
            xbeta = dN1dbeta*x1 + dN2dbeta*x2 + dN3dbeta*x3 + dN4dbeta*x4 + dN5dbeta*x5 + dN6dbeta*x6 + dN7dbeta*x7 + dN8dbeta*x8
            xgamma = dN1dgamma*x1 + dN2dgamma*x2 + dN3dgamma*x3 + dN4dgamma*x4 + dN5dgamma*x5 + dN6dgamma*x6 + dN7dgamma*x7 + dN8dgamma*x8
            yalpha = dN1dalpha*y1 + dN2dalpha*y2 + dN3dalpha*y3 + dN4dalpha*y4 + dN5dalpha*y5 + dN6dalpha*y6 + dN7dalpha*y7 + dN8dalpha*y8
            ybeta = dN1dbeta*y1 + dN2dbeta*y2 + dN3dbeta*y3 + dN4dbeta*y4 + dN5dbeta*y5 + dN6dbeta*y6 + dN7dbeta*y7 + dN8dbeta*y8
            ygamma = dN1dgamma*y1 + dN2dgamma*y2 + dN3dgamma*y3 + dN4dgamma*y4 + dN5dgamma*y5 + dN6dgamma*y6 + dN7dgamma*y7 + dN8dgamma*y8
            zalpha = dN1dalpha*z1 + dN2dalpha*z2 + dN3dalpha*z3 + dN4dalpha*z4 + dN5dalpha*z5 + dN6dalpha*z6 + dN7dalpha*z7 + dN8dalpha*z8
            zbeta = dN1dbeta*z1 + dN2dbeta*z2 + dN3dbeta*z3 + dN4dbeta*z4 + dN5dbeta*z5 + dN6dbeta*z6 + dN7dbeta*z7 + dN8dbeta*z8
            zgamma = dN1dgamma*z1 + dN2dgamma*z2 + dN3dgamma*z3 + dN4dgamma*z4 + dN5dgamma*z5 + dN6dgamma*z6 + dN7dgamma*z7 + dN8dgamma*z8
            
            coordlist = [xalpha, xbeta, xgamma, yalpha, ybeta, ygamma, zalpha, zbeta, zgamma]

            # Zero out values close to zero
            threshold = 1e-10

            zeroout = []
            for a in coordlist:
                if a < threshold:
                    zeroout.append(0.0)
                else:
                    zeroout.append(a)

            # Natural coordinates
            xs, xt, xu, ys, yt, yu, zs, zt, zu = zeroout

            # Jacobian - maps natural coordinates to physical 
            J = np.array([[xs, xt, xu], 
                    [ys, yt, yu],
                    [zs, zt, zu]])

            detJ = xalpha*(ybeta*zgamma - zbeta*ygamma) - yalpha*(xbeta*zgamma - zbeta*xgamma) + zalpha*(xbeta*ygamma - ybeta*xgamma)

            rightmatrix = np.array([[dN1dalpha, dN2dalpha, dN3dalpha, dN4dalpha, dN5dalpha, dN6dalpha, dN7dalpha, dN8dalpha], 
                                    [dN1dbeta, dN2dbeta, dN3dbeta, dN4dbeta, dN5dbeta, dN6dbeta, dN7dbeta, dN8dbeta], 
                                    [dN1dgamma, dN2dgamma, dN3dgamma, dN4dgamma, dN5dgamma, dN6dgamma, dN7dgamma, dN8dgamma]])

            # Shape function derivates 
            Nxyz = np.linalg.solve(J, rightmatrix)

            N1x, N2x, N3x, N4x, N5x, N6x, N7x, N8x = Nxyz[0][0:8]
            N1y, N2y, N3y, N4y, N5y, N6y, N7y, N8y = Nxyz[1][0:8]
            N1z, N2z, N3z, N4z, N5z, N6z, N7z, N8z = Nxyz[2][0:8]

            # Extrapolated strain matrix
            B = np.array([[N1x, 0.0, 0.0, N2x, 0.0, 0.0, N3x, 0.0, 0.0, N4x, 0.0, 0.0, N5x, 0.0, 0.0, N6x, 0.0, 0.0, N7x, 0.0, 0.0, N8x, 0.0, 0.0], [0.0, N1y, 0.0, 0.0, N2y, 0.0, 0.0, N3y, 0.0, 0.0, N4y, 0.0, 0.0, N5y, 0.0, 0.0, N6y, 0.0, 0.0, N7y, 0.0, 0.0, N8y, 0.0], [0.0, 0.0, N1z, 0.0, 0.0, N2z, 0.0, 0.0, N3z, 0.0, 0.0, N4z, 0.0, 0.0, N5z, 0.0, 0.0, N6z, 0.0, 0.0, N7z, 0.0, 0.0, N8z], [N1y, N1x, 0.0, N2y, N2x, 0.0, N3y, N3x, 0.0, N4y, N4x, 0.0, N5y, N5x, 0.0, N6y, N6x, 0.0, N7y, N7x, 0.0, N8y, N8x, 0.0], [0.0, N1z, N1y, 0.0, N2z, N2y, 0.0, N3z, N3y, 0.0, N4z, N4y, 0.0, N5z, N5y, 0.0, N6z, N6y, 0.0, N7z, N7y, 0.0, N8z, N8y], [N1z, 0.0, N1x, N2z, 0.0, N2x, N3z, 0.0, N3x, N4z, 0.0, N4x, N5z, 0.0, N5x, N6z, 0.0, N6x, N7z, 0.0, N7x, N8z, 0.0, N8x]])
        
            BE_D = np.dot(np.transpose(B), D)
            BEB_D = np.dot(BE_D, B) * detJ

            # Update stiffness matrix
            ke = ke + BEB_D

            return ke

    def getElementStiffness(self):
        
        Ae1 = np.zeros((self.topology.getNumElements(), self.topology.getNumElements()))
        Ke1 = np.zeros((24, 24, self.topology.getNumElements()))

        for element in range(self.topology.getNumElements()):
            elementcoord = self.coord[self.con[element].conj().transpose()]
            Ae1[element, 0] = ConvexHull(elementcoord).vertices[0][0]
            Ke1[:, :, element] = self.hexElementStiffness(element)
             
        return Ke1, Ae1

    def calcEDOFS(self):
        
        for element in range(self.topology.getNumElements()):

            self.edofs.append(np.array([3*self.con[element][0]-3, 3*self.con[element][0]-2, 3*self.con[element][0]-1,
                3*self.con[element][1]-3, 3*self.con[element][1]-2, 3*self.con[element][1]-1,
                3*self.con[element][2]-3, 3*self.con[element][2]-2, 3*self.con[element][2]-1,
                3*self.con[element][3]-3, 3*self.con[element][3]-2, 3*self.con[element][3]-1,
                3*self.con[element][4]-3, 3*self.con[element][4]-2, 3*self.con[element][4]-1,
                3*self.con[element][5]-3, 3*self.con[element][5]-2, 3*self.con[element][5]-1,
                3*self.con[element][6]-3, 3*self.con[element][6]-2, 3*self.con[element][6]-1,
                3*self.con[element][7]-3, 3*self.con[element][7]-2, 3*self.con[element][7]-1]).astype(int))

    def getEDOFS(self):
        return self.edofs


    def stiffness(self, Ke, design):

        K = sp.lil_matrix((3*self.topology.getNumConnections(), 3*self.topology.getNumConnections()))
        for element in range(self.topology.getNumElements()):
            K[np.ix_(self.edofs[element], self.edofs[element])] += Ke[:, :, element] * design[element]

        return K.tocsc()      