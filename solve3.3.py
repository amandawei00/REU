import sys
import numpy as np
import scipy.integrate as integrate
import csv
import matplotlib.pyplot as plt

# sys.path.append("DSS14_Python")
sys.path.append("python_scripts")
sys.path.append("DSSLIB")

# from DSS14_Python import DSS
from DSS_Python import DSS
import lhapdf as pdf
from bk_interpolate import N

# MAKE USER FRIENDLY
# - IH FOR HADRON TYPE, IC FOR HADRON CHARGE SHOULD BE MODIFIABLE AND INITIATED UPON CONSTRUCTION

class Master():
    def __init__(self,y,s_NN,qsq,K,h):
        self.p = pdf.mkPDF("CT10",0)
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

        qual = ''
        if self.f == -3: 
            i = 5
            qual = 'anti-strange quark'
        elif self.f == -2: 
            i = 1
            qual = 'anti-up quark'
        elif self.f == -1: 
            i = 3
            qual = 'anti-down quark'
        elif self.f == 1: 
            i = 2
            qual = 'down quark'
        elif self.f == 2: 
            i = 0
            qual = 'up quark'
        elif self.f == 3: 
            i = 4
            qual = 'strange'
        else:
            return 0.0

        pdf_qp = self.p.xfxQ2(self.f,x1,q2) # returns x1*f(x1,pt^2) where f is pdf
        bkf = self.n.udg_f(x2,q/z)
        ff_hq = self.ff.get_f(z,q2,self.hadron)[i]

        # print("pdf_q = " + str(pdf_qp) + ", bkf = " + str(bkf) + ", ff_hq = " + str(ff_hq))
        print("ff_hq = " + str(ff_hq) + ", z = " + str(z) + ", q2 = " + str(q2))
        a = (1/np.power(z,2))*(pdf_qp*bkf*ff_hq)

        return a
###################################################################################################################################
    def integrand1(self, z):
        xf = self.get_xf()

        x1 = xf/z
        x2 = x1*np.exp(-2*self.yh)
        q = self.p_t
        q2 = q*q

        # if self.f == 21
        i=8 # for gluons

        pdf_gp = self.p.xfxQ2(self.f,x1,q2)
        bka = self.n.udg_a(x2,q/z)
        ff_hg = self.ff.get_f(z,q2,self.hadron)[i]
        return (1/np.power(z,2))*(pdf_gp*bka*ff_hg)

################################################################################################################################
    def rhs(self,pt): # DEFINE TOL
        self.p_t = pt
        xf = self.get_xf()
        m = 0.0

        gluon_intg = integrate.quad(self.integrand1, xf, 1.0)[0]
        for i in self.flavors:
            self.f = i
            quark = integrate.quad(self.integrand,xf,1.0)[0]
            print("quark = " + str(quark) + ", gluon = " + str(gluon_intg))
            m += quark + gluon_intg # integral

        return m*self.K/np.power(2*np.pi, 2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if __name__=="__main__":
    ih = 'pi0' # hadron type
    y = 3.3
    s_NN = np.power(200,2) # GeV
    qsq2 = 0.4
    K = 1.0

    s = Master(y, s_NN, qsq2, K, ih)
   
    n = 5
    a = 1.01
    b = 5.0
    dp_t = (b - a)/n

    p_t = np.arange(a,b,dp_t)    
    cs = np.zeros(len(p_t))
    for i in range(len(p_t)):
        cs[i] = s.rhs(p_t[i])
        print(str(p_t[i])+", "+str(cs[i]))


    with open('output_pi0_y-3.3_x0-0.02.csv', "w") as csvfile:
        writer = csv.writer(csvfile, delimiter = '\t')
        
        for i in range(len(p_t)):
            writer.writerow([p_t[i],cs[i]])

    plt.plot(p_t, cs)
    plt.show()
