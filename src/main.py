from optimize import Optimizer
from topology import Topology

if __name__ == "__main__":
    
    nelx, nely, nelz = 40, 2, 2
    design = Topology(nelx, nely, nelz, 40, 10, 20)

    forceLocation = 3 * ((nelx)*(nely+1)*(nelz+1)+(nely/2+1)*(nelz+1))

    optimizer = Optimizer(design, 0.15, -1000, forceLocation, 200e9)