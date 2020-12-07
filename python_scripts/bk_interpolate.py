# from FBT import FBT
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import scipy.integrate as intg
import time


class N:

    def __init__(self):
        self.n_ = 400
        self.x0 = 0.02
        self.xr1 = np.log(3.e-6)
        # self.xr1 = 0.0
        self.xr2 = np.log(60.e0)  # limit of integration in fourier transform calculation
        tol = 1.e-8
        self.width = 100


        # read results.csv file from BK solution to pandas dataframe
        self.df = pd.read_csv("full_bk_results.csv", sep="\t")
        self.df.columns = ['kuta', 'y', 'vr', 'vfr', 'prev']
        # self.df = self.df.drop(self.df[self.df.kuta == "kuta"].index)

        # converting dataframe element types
        self.df["kuta"] = self.df["kuta"].astype('int')
        self.df["y"] = (self.df["y"].astype('float32')).round(decimals=1)
        self.df["vr"] = self.df["vr"].astype('float64')
        self.df["vfr"] = self.df["vfr"].astype('float64')

        self.df = self.df.loc[(self.df['kuta'] == 4)]  # only consider 4th-order RK solutions

        # self.r = np.unique(self.df.vr.to_numpy(dtype='float64'))
        # nump_vr = np.array(self.df.vr, dtype=np.float64)
        # self.r = list(dict.fromkeys(nump_vr))
        self.r = np.concatenate(np.array(self.df.loc[self.df['y'] == 0.][['vr']]))  # r values over which N(r,Y) are evaluated
        self.y_vals = np.unique(np.array(self.df[['y']]))

        self.interp_y = []  # to be filled by interpolated functions over r at some fixed y
        for i in self.y_vals:
            sub_ = self.df.loc[(self.df['kuta'] == 4) & (self.df['y'] == i.round(decimals=1))]
            x_ = np.concatenate(np.array(sub_[['vr']]), axis=0)
            y_ = np.concatenate(np.array(sub_[['vfr']]), axis=0)

            self.interp_y.append(interp1d(x_, y_, kind='cubic'))

        self.interp_r = []  # to be filled by interpolated functions over y at some fixed r
        for r_ in self.r:
            sub_ = self.df[(self.df['vr'] < (r_ + tol)) & (self.df['vr'] > (r_ - tol))]  # all_y = self.df['y']
            x_ = np.concatenate(np.array(sub_[['y']]), axis=0)  # np.concatenate is needed because np.array returns an array of arrays
            y_ = np.concatenate(np.array(sub_[['vfr']]), axis=0)

            # print(x_)
            # print(y_)
            # x_ = np.concatenate(sub_[['y']].to_numpy(), axis=0)
            # y_ = np.concatenate(sub_[['vfr']].to_numpy(), axis=0)
            self.interp_r.append(interp1d(x_, y_, kind='cubic'))

    def bk_f(self, y):  # fundamental representation
        t1 = time.time()
        rap = np.around(np.arange(-4.9, 30.0, 0.1), 2)  # must round np.arange values to 2 decimal places
        # search file to see if interpolation already exists !!!
        if y in rap:
            i = np.where(rap == y)[0][0]
            alpha = self.interp_y[i]
        else:
            index = np.where(rap == round(y, 1))
            a = np.zeros(len(self.r))  # to be filled with N(r,y) values over which to interpolate
            for i in range(self.width):
                f = self.interp_r[i]
                a[i] = f(y)
            alpha = interp1d(self.r, a, kind='cubic')


        t2 = time.time()
        print("bk_f took " + str(t2-t1) + " seconds")
        return alpha




        # write interpolated object to file !!!

    def bk_a(self, y):
        t1 = time.time()
        f = self.bk_f(y)
        x = self.r

        g = 2 * f(x) - np.power(f(x), 2)
        h = interp1d(self.r, g, kind='cubic')

        t2 = time.time()
        print("bk_a took " + str(t2-t2) + " seconds")
        return h

    def udg_f(self, x, k):
        t1 = time.time()
        f = self.bk_f(np.log(self.x0 / x))
        integrand = lambda r: (1 - f(r)) * self.bessel(k * r, 0)
        a = 2 * np.pi * intg.quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=1.e-4)[0]
        t2 = time.time()
        print("udg_f took " + str(t2-t1) + " seconds")
        return a

    def udg_a(self, x, k):
        t1 = time.time()
        f = self.bk_a(np.log(self.x0 / x))
        integrand = lambda r_: (1 - f(r_)) * self.bessel(k * r_, 0)
        a = 2 * np.pi * intg.quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=1.e-4)[0]
        t2 = time.time()
        print("udg_a took " + str(t2-t1) + " seconds")
        return a

    def bessel(self, x, alpha):
        # print("entering bessel")
        f = lambda t: np.cos(alpha * t - x * np.sin(t))
        g = (1 / np.pi) * intg.quad(f, 0, np.pi, epsabs=1.e-3)[0]
        # print("exiting bessel")
        return g

# end of class


if __name__ == "__main__":
    n = N()

    print(n.udg_a(0.001161195799828, 2.0959615528))

