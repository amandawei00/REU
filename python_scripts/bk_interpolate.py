import numpy as np
import scipy.interpolate as interpolate
import scipy.integrate as intg
import pandas as pd
import csv

class N():
    def __init__(self):
        self.n_ = 400
        self.x0 = 0.02
        self.xr1 = np.log(3.e-6)
        self.xr2 = np.log(60.e0) # limit of integration in fourier transform calculation

        # read results.csv file from BK solution to pandas dataframe
        self.df = pd.read_csv("results.csv", sep="\s+")
        self.df = self.df.drop(self.df[self.df.kuta == "kuta"].index)

        # converting dataframe element types
        self.df["kuta"] = self.df["kuta"].astype('int')
        self.df["y"] = (self.df["y"].astype('float32')).round(decimals=1)
        self.df["vr"] = self.df["vr"].astype('float64')
        self.df["vfr"] = self.df["vfr"].astype('float64')

        self.r = self.df.vr.to_numpy(dtype = 'float64')
        self.y = self.df.y.to_numpy(dtype = 'float32')
        self.n = self.df.vfr.to_numpy(dtype = 'float64')

        self.funcs = []
        for i in np.arange(0.0,30.0,0.10):
            sub_ = self.df.loc[(self.df['kuta'] == 4.0) & (self.df['y'] == i)]

            x = sub_[['vr']].to_numpy()

    def bk_f(self,y): # fundamental representation
        # search file to see if interpolation already exists !!!

        # isolating part of df containing r, N(r,Y) values for rapidity scale of interest
        sub_ = self.df.loc[(self.df['kuta'] == 4.0) & (self.df['y'] == y)]

        # x_ = np.concatenate(sub_[['vr']].to_numpy(), axis = 0)
        # y_ = np.concatenate(sub_[['vfr']].to_numpy(), axis = 0)
    
        x_ = sub_[['vr']].to_numpy()
        y_ = sub_[['vfr']].to_numpy()

        f = interpolate.CubicSpline(x_, y_) # creates interpolation

        # write interpolated object to file !!!

        return f

    def bk_a(self, y):
        f = self.bk_f(y)
        return 2*f - np.power(f,2)

    def udg_f(self, x,  k):
        f = self.bk_f(np.log(self.x0/x))
        integrand = lambda r: (1 - f(r))*self.bessel(k*r,0)
        return 2*intg.quad(integrand, 0, self.xr2)[0]

    def udg_a(self, x, k):
        f = self.bk_a(np.log(self.x0/x))
        integrand = lambda r: (1-f(r))*self.bessel(k*r,0)
        return 2*intg.quad(integrand, 0, self.xr2)[0]

    def bessel(self, x, alpha):
        f = lambda t: np.cos(alpha*t - x*np.sin(t))
        return (1/np.pi)*intg.quad(f,0,np.pi)[0]

    # def bessel_intg(self,k):
    #    f = lambda u: (u/np.power(k,2)) * self.bessel(u,0)
    #    return intg.quad(f,0.0,np.inf)

# end of class

if __name__ == "__main__":
    n = N()
