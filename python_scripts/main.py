import numpy as np
import scipy.interpolate as interpolate
import scipy.integrate as integrate

class Test():
    def __init__(self):
        self.n=399
        self.xr1=np.log(3.e-6)
        self.xr2=np.log(60.e0)
        self.hy=0.1
        self.ymax=30.0
        self.y=np.arange(0.0,self.ymax+1,self.hy)

        self.r = np.zeros(400)
        self.xlr = np.zeros(400)
        self.n_current = np.zeros(400)

        self.xprec1=5e-3
        self.xprec2=5e-4
        self.xprec3=1.5e-2

        self.nmax1=2000
        self.nmax2=20000
        self.rlim1=5e-8
        self.rlim2=0.9999

        self.xlam=0.241
        self.af=1.0
        self.ge=0.5772156659
        self.xnf=3.0
        self.beta=(33-2.0*self.xnf)/(12.0*np.pi)
        self.pre=0.5/self.beta
        self.auxal=2.0/self.xlam
        self.iglob=1

        self.xn=60.e0
        self.qsq2=0.168e0*self.xn
        self.gamma=1.0

        self.method="RK$"

        for i in range(400):
            self.xlr[i]=self.xr1+i*self.hy
            self.r[i]=np.exp(self.xlr[i])
            xlog=np.log(1.0/(0.241*self.r[i])+np.exp(1.0))
            xexp=np.power(self.qsq2*np.power(self.r[i],2),self.gamma)*xlog/4.0
            self.n_current[i] = 1.0-np.exp(-xexp)

    def r2(self,r,r1,theta):
        r2=np.sqrt(np.power(r,2)+np.power(r1,2)-2*r*r1*np.cos(theta))
        return r2

# atm, consider running coupling frozen at alpha=0.7 UNFINISHED
    def alphaS(self,r):
        return 0.7

    def kernel(self,rvar,r1,theta):
        r2=self.r2(rvar,r1,theta)

        Nc=3.0
        factor=Nc*self.alphaS(rvar)/(2*np.power(np.pi,2))
        term1=np.power(rvar,2)/(np.power(r1,2)*np.power(r2,2))
        term2=(1/np.power(r1,2))*(self.alphaS(r1)/self.alphaS(r2)-1)
        term3=(1/np.power(r2,2))*(self.alphaS(r2)/self.alphaS(r1)-1)

        return factor*(term1+term2+term3)

    def n(self,f5,x):
        n_=0.0
        if x<self.xr1:
            ex=1.99
            n_=np.power(self.n_current[0]*np.exp(x)/self.r[0],ex)
        elif x>self.xr2:
            n_=1.0
        else:
#            f5=interpolate.interp1d(self.xlr,self.n_current,kind=5)
            n_=f5(x)

        if n_<0.0: return 0.0
        if n_>1.0: return 1.0

        return n_

# integrate function f for inner bounds [a,b], outer bounds [c,d], to precision prec
    # check to see if this works
    def integrate(self, f, a, b,  tol):
        # weights and abscissas for Gaussian quadrature at n=8:
        w = [0.101228536290376259, 0.222381034453374471, 0.313706645877887287, 0.362683783378361983,
             0.0271524594117540949, 0.0622535239386478929, 0.0951585116824927848, 0.124628971255533872,
             0.149595988816576732, 0.169156519395002538, 0.182603415044923589, 0.189450610455068496]
        x = [0.960289856497536232, 0.796666477413626740, 0.525532409916328986, 0.183434642495649805,
             0.989400934991649933, 0.944575023073232576, 0.865631202387831744, 0.755404408355003034,
             0.617876244402643748, 0.458016777657227386, 0.281603550779258913, 0.0950125098376374402]

        s = 0.0
        iter=0.0
        cont1=True
        if b == a: return s
        const = 0.005 / (b - a)
        bb = a
        iter += 1
        while cont1:
            aa = bb
            bb = b

            cont = True
            iter = 0
            while cont:
                c1 = 0.5 * (bb + aa)
                c2 = 0.5 * (bb - aa)
                s8 = 0.0

                for i in range(5):
                    u = c2 * x[i]
                    s8 += w[i] * (f(c1 + u) + f(c1 - u))

                s8 = c2 * s8
                s16 = 0.0

                for i in range(5, 12):
                    u = c2 * x[i]
                    s16 += w[i] * (f(c1 + u) + f(c1 - u))

                s16 = c2 * s16

                if np.abs(s16 - s8) < tol * np.abs(s16): cont = False
                bb = c1
                if 1.0 + np.abs(const * c2) != 1.0: cont = True

            s += s16
            if iter > self.nmax2: iter = 0
            if bb != b:
                cont1 = True
            else:
                cont1 = False
        return s

    def inner(self,f5,r,r1,t):
        r1 = np.exp(r1)
        q1 = 0.25*np.power(r,2)+np.power(r1,2)-r*r1*np.cos(t)
        q2 = 0.25*np.power(r,2)+np.power(r1,2)+r*r1*np.cos(t)

        if q1 < 0.0: q1=0.0
        if q2 < 0.0: q2=0.0

        fr = f5(r) #!!!!!!!!!!!!!!CHECK!!!!!!!!!!!!!#
        fr1 = self.n(f5,np.log(np.sqrt(q1)))
        fr2 = self.n(f5,np.log(np.sqrt(q2)))

        ker = self.kernel(r,r1,t)

        return np.power(r1,2)*ker*(fr1+fr2-fr-fr1*fr2)


    def xk(self,f5):
        integral=np.zeros(400)
        xp1=5.e-3
        xp2=5.e-3

        for i in range(len(self.r)):
            integral[i]=self.integrate(lambda t_: self.integrate(lambda r1_: self.inner(f5,self.r[i],r1_,t_),self.xr1,self.xr2,xp1),0,0.5*np.pi,xp2)

        for i in range(self.n):
            print(integral[i])
        return integral

    def evolve(self):
        # 4th order runge-kutta, and then call interpolations after correction for self.n
        func=lambda x: interpolate.interp1d(self.xlr,self.n_current,kind=5)
        for i in range(len(self.y)):
            y0=self.y[i]
            print("Y=",y0)

            if self.method == "Euler": kkuta = 1
            elif self.method == "RK2": kkuta = 2
            elif self.method == "RK4": kkuta = 4

            n_correction=self.xk(func)
            n_corrected=self.n_current

            dy = 0.1
            for k in np.arange(1,kkuta+1):
                if kkuta == 1: dy = self.hy
                elif kkuta == 2:
                    if k == 1: dy = 0.5*self.hy
                    if k == 2: dy = self.hy
                elif kkuta == 4:
                    if k == 1: dy = 0.5*self.hy
                    if k == 2: dy = 0.5*self.hy
                    if k == 3: dy = self.hy
                    if k == 4: dy = self.hy

                for j in range(len(self.r)):
                    n_corrected[j] += dy*n_correction[j]
                    if n_corrected[j] > 1.0: n_corrected[j] = 1.0

                n_next=np.zeros(self.n)
                for i in range(len(self.r)):
                    n_next[i] = self.n_current[i]+self.hy*n_correction[i]
                    if n_next[i] > 1.0: n_next[i] = 1.0

                self.n_current = n_next


if __name__ == "__main__":
    t = Test()
    t.evolve()




