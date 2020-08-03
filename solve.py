import sys
import numpy as np
sys.path.append("DSS14_Python")
sys.path.append("python_scripts")

from DSS14 import DSS
import lhapdf as pdf
from bk_interpolate import N

# MAKE USER FRIENDLY
# - IH FOR HADRON TYPE, IC FOR HADRON CHARGE SHOULD BE MODIFIABLE AND INITIATED UPON CONSTRUCTION
# - FIX FRAGMENTATION FUNCTION
# - FINISH CODE FOR BK INTERPOLATION

class Master():
    def __init__(self,y,s_NN,qsq,K,h):
        self.p = pdf.makePDF("CT10",0)
        self.p = pdf.mkPDF("CT10/0")

        self.n = N()
        self.ff = DSS()

        self.qsq2 = qsq # saturation scale
        self.sNN = s_NN # collision energy per nucleon [GeV]
        self.K = K
        self.hadron = h

        # self.flavors = [1,2,3,4,5,21] # flavors: d(1), u(2), s(3), c(4), b(5), g(21)
        self.flavors = self.p.flavors()
        self.f = 0.0

        self.p_t = 0.0 # plotting points w respect to p_t
        self.yh = y
 ###################################################################################################################################
    def get_xf(self):
        xf = (self.p_t/np.sqrt(self.sNN))*np.exp(self.yh)
        return xf
####################################################################################################################################
    def integrand(self,z):
        xf = self.get_xf()

        x1 = xf/z
        x2 = x1*np.exp(-2*self.yh)
        q = self.p_t
        q2 = q*q

        if self.f == -3: i = 5
        elif self.f == -2: i = 1
        elif self.f == -1: i = 3
        elif self.f == 1: i = 2
        elif self.f == 2: i = 0
        elif self.f == 3: i = 4
        elif self.f == 21: i = 8
        else:
            return 0.0

        pdf_qp = self.p.xfxQ2(self.f,x1,q2) # returns x1*f(x1,pt^2) where f is pdf
        bkf = self.n.udg_f(x2,q/z)
        # ff_hq = self.ff.fDSS(4,-1,0,z,q2)[i] # for negatively charged hadron DEFINE i HERE SPECIFIC TO COLOR
        ff_hq = self.ff.get_f(z,q2,self.hadron)[i]

        pdf_gp = self.p.xfxQ2(21,x1,q2) # 21 for gluon
        bka = self.n.udg_a(x2,q/z)
        # ff_hg = self.ff.fDSS(4,-1,0,z,q2)[8]
        ff_hg = self.ff.get_f(z,q2,self.hadron)[8]

        a = (1/np.power(z,2))*(pdf_qp*bkf*ff_hq + pdf_gp*bka*ff_hg)

        return a
###################################################################################################################################
    def rhs(self,pt): # DEFINE TOL
        self.p_t = pt
        xf = self.get_xf()
        m = 0.0

        for i in self.flavors():
            self.color = i
            m += integrate.quadrature(self.integrand, xf, 1.0, tol) # integral

        return m*self.K/np.power(2*np.pi, 2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if __name__=="__main__":
   ih = 'pi0' # hadron type
   y = 3.3
   s_NN = np.power(200,2) # GeV
   qsq2 = 0.4
   K = 0.4

   s = Master(y, s_NN, qsq2, K, ih)
   print(s.rhs(2.0))
