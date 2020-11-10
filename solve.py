import sys
import numpy as np
import scipy.integrate as integrate
import csv
import matplotlib.pyplot as plt
import subprocess

# sys.path.append("DSS14_Python")
sys.path.append("python_scripts")
# sys.path.append("DSSLIB")

# from DSS14_Python import DSS14
# from DSS_Python import DSS
import lhapdf as pdf
from bk_interpolate import N

# - IH FOR HADRON TYPE, IC FOR HADRON CHARGE SHOULD BE MODIFIABLE AND INITIATED UPON CONSTRUCTION

class Master():
    def __init__(self, h, y, qsq, s_NN, K):
        self.p = pdf.mkPDF("CT10",0)
        self.p = pdf.mkPDF("CT10/0")
# 	self.p = pdf.mkPDFs("CT10nlo",0)

        self.n = N()
        self.ff = pdf.mkPDF("DSS07PI",0)

	print("done done ")
        self.qsq2 = qsq # saturation scale
        self.sNN = s_NN # collision energy per nucleon [GeV]
        self.K = K
        self.hadron = h

        # self.flavors = [1,2,3,4,5,21] # flavors: d(1), u(2), s(3), c(4), b(5), g(21)
        self.flavors = self.p.flavors()
        self.f = 0.0

        self.p_t = 0.0 # plotting points w respect to p_t
        self.yh = y
	
	print("class successfully instantiated")
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

        """qual = ''
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
            return 0.0 """

        pdf_qp = self.p.xfxQ2(self.f,x1,q2) # returns x1*f(x1,pt^2) where f is pdf
        bkf = self.n.udg_f(x2,q/z)
        ff_hq = self.ff.xfxQ2(self.f, z, q2)

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
        ff_hg = self.ff.xfxQ2(self.f, z,q2)
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
            # print("quark = " + str(quark) + ", gluon = " + str(gluon_intg))
            m += quark + gluon_intg # integral

        return m*self.K/np.power(2*np.pi, 2)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if __name__=="__main__":

    param = []
    """ order of parameters in parameters.txt:
		1. hadron type
		2. y
		3. qsq2
		4. s_NN
		5. K
		6. p_t (lower bound)
		7. p_t (upper bound)
		8. number of points evaluated """

    with open("parameters.txt", "r") as infile:  # opening parameters file to read in input
	reader = csv.reader(infile)
	
	for r in reader:  
	    raw = r[0].split()
	    param.append(raw[len(raw)-1])  # reading in values

    s = Master(param[0], float(param[1]), float(param[2]), float(param[3]), float(param[4]))  # creating instance of class
   
    a = float(param[5])
    b = float(param[6])
    n = int(param[7])
    dp_t = (b - a)/n

    p_t = np.arange(a,b,dp_t)  # values of tranverse momenta over which differential cross section will be evaluated
    cs = np.zeros(len(p_t))
    
    with open("temp.csv", "w") as tfile:  # write temporary output file so progress can be checked without interrupting code
	writer = csv.writer(tfile, delimiter='\t')

    	for i in range(len(p_t)):
            cs[i] = s.rhs(p_t[i])
            writer.writerow([p_t[i],cs[i]])


    fname = "run_1"  # folder name
    subprocess.run(["mkdir", "output/" + fname])
    subprocess.run(["cp", "parameters.txt", "output/" + fname])  # copy parameter file to output folder
    subprocess.run(["cd", "output/" + fname])  # enter directory to write output here
	
    with open("results.csv", "w") as csvfile:
        writer = csv.writer(csvfile, delimiter = '\t')
        
        for i in range(len(p_t)):
            writer.writerow([p_t[i],cs[i]])

    subprocess.run(["cd", "../.."])  # return to original directory
    
    plt.plot(p_t, cs)
    plt.show()

# end of program
