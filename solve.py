import sys
import numpy as np
sys.path.append("DSS14_Python")
sys.path.append("python_scripts")

from DSS14 import DSS as ff
import lhapdf as pdf
from bk_interpolate import N


class main():
    def __init__(self,y,s_NN,fl):
        self.p = pdf.makePDF("CT10",0)
        self.p = pdf.mkPDF("CT10/0")

        self.q = 0.4 # saturation scale
        self.sNN = s_NN # collision energy per nucleon [GeV]
        self.yh = y

        self.flavors = [1,2,3,4,5,21] # flavors: d(1), u(2), s(3), c(4), b(5), g(21)
        self.q = 0.0
        self.n = N() # to find N, run n.bkf(r,Y)
        
    
    def integrand(self,z):
        x1 = self.xf/z
        x2 = x1*np.exp(-2*self.yh)

        # t1 = x1*self.p.xfxQ2(self.q,x1,np.power(pt,2))

    def integral(self,pt):
        xF = pt*np.exp(self.yh)/np.sqrt(self.sNN)

        # return integrate(integrand,xF,1,prec)


        
if __name__=="__main__":
    s_NN = 200 # GeV
    y_h1 = 2.2 # pseudorapidities for negative-charge hadrons
    y_h2 = 3.2 

    y_h3 = 3.3 #pseudorapidities for neutral pions
    y_h4 = 3.8 
    y_h5 = 4.0
    
    qsq2 = 0.4 # GeV^2 quark saturation scale for gold nucleus
    x0 = 0.02 
