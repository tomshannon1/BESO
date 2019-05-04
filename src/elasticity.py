import numpy as np

class Elasticity:

    def __init__(self, ym, pr):
        
        self.ym = ym
        self.pr = pr
        
    def getElasticityMatrix(self):

        scalar = (self.ym / ((1+self.pr)*(1-2*self.pr)))        
        return scalar * np.matrix([[1-self.pr, self.pr, self.pr, 0, 0, 0],
                          [self.pr, 1-self.pr, self.pr, 0, 0, 0],
                          [self.pr, self.pr, 1-self.pr, 0, 0, 0],
                          [0, 0, 0, (1-2*self.pr)/2, 0, 0],
                          [0, 0, 0, 0, (1-2*self.pr)/2, 0],
                          [0, 0, 0, 0, 0, (1-2*self.pr)/2]])
                          