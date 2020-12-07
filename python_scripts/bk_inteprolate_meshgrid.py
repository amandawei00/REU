# from FBT import FBT
import numpy as np
import pandas as pd
import scipy.interpolate as interp
import scipy.integrate as intg
import time
import cmath


class N:

    def __init__(self):
        self.n_ = 400
        self.x0 = 0.02
        self.xr1 = np.log(3.e-6)
        # self.xr1 = 0.0
        self.xr2 = np.log(60.e0)  # limit of integration in fourier transform calculation
        tol = 1.e-8
        self.width = 100
        self.pointsy = 400
        self.pointsr = 600


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

        self.r = np.concatenate(np.array(self.df.loc[self.df['y'] == 0.][['vr']]))  # r values over which N(r,Y) are evaluated
        self.y = np.unique(np.array(self.df[['y']]))

        points = [(jingle, bell) for jingle in self.r for bell in self.y]  # coordinates made from permutations of values in self.r, self.y
        values = [self.df.loc[(self.df['y'] == p[0]) & ((self.df['vr'] < (p[1]+tol)) & (self.df['vr'] > (p[1]-tol)))][['vfr']].iloc[0]['vfr'] for p in points]  # values of N(r,Y) which we KNOW from BK solution
        self.ri = np.mgrid[self.r[0]:self.r[len(self.r)-1]:complex(self.pointsr)]  # grid over r
        self.yi = np.mgrid[self.y[0]:self.y[len(self.y)-1]:complex(self.pointsy)]  # grid over y
        self.f = self.interp.griddata(points, values, (self.ri, self.yi), method='cubic')  # interpolated grid
        print("grid created, interpolated")

    # finds location of val in grid in the sense that if val is between two elements in grid, find_index will return the index of the lower element
    def find_index(self, val, grid):
        index = 0

        cont = True
        i = 0
        while cont:
            if val >= grid[i]:
                index = i
                cont = False

            i += 1

            if i >= len(grid):
                cont = False

        return index

    # given any r, Y within bounds of self.r and self.y, master returns interpolated value of N(r,Y)
    def master(self, r, y):
        rl = self.find_index(r, self.ri)  # LOCATION of r in self.ri
        yl = self.find_index(y, self.yi)  # LOCATION of y in self.yi

        n = 0.0

        if (r == self.ri[len(self.ri)-1]) and (y == self.yi[len(self.yi)-1]):
            n = self.f[rl][yl]
        elif rl >= len(self.ri)-1:  # edge case: do 1d interpolation over y
            xvals = [self.yi[yl], self.yi[yl+1]]
            yvals = [self.f[rl][yl], self.f[rl][yl+1]]
            n = interp.interp1d(xvals, yvals, kind = 'linear')

        elif yl >= len(self.yi)-1:  # edge case: do 1d interpolation over r
            xvals = [self.ri[rl], self.ri[rl+1]]
            yvals = [self.f[rl][yl], self.f[rl+1][yl]]
            n = interp.interp1d(xvals, yvals, kind = 'linear')
        else:
            xvals = [self.ri[rl], self.ri[rl+1]]
            yvals = [self.yi[yl], self.yi[yl+1]]
            zvals = [[self.f[rl][yl], self.f[rl+1][yl]],[self.f[rl][yl+1], self.f[rl+1][yl+1]]]

            n = interp.interp2d(xvals, yvals, zvals, kind='linear')

        return n

    def master_adj(self, r, y):
        return 2 * self.master(r, y) - np.power(self.master(r, y), 2)

    def udg_f(self, x, k):
        t1 = time.time()
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.master(r_, y_)) * self.bessel(k * r_, 0)
        a = 2 * np.pi * intg.quad(integrand, self.r[0], self.r[len(self.r)-1], epsabs=1.e-4)[0]
        t2 = time.time()
        print("udg_f took " + str(t2-t1) + " seconds")
        return a

    def udg_a(self, x, k):
        t1 = time.time()
        y_ = np.log(self.x0 / x)
        integrand = lambda r_: (1 - self.master_adj(r_, y_)) * self.bessel(k * r_, 0)
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