# sys.path.append("Ogata/python")
# from FBT import FBT
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
import scipy.integrate as intg

class N():
    def __init__(self):
        self.n_ = 400
        self.x0 = 0.02
        # self.xr1 = np.log(3.e-6)
        self.xr1 = 0.0
        self.xr2 = np.log(60.e0)  # limit of integration in fourier transform calculation

        # read results.csv file from BK solution to pandas dataframe
	self.df = pd.read_csv("full_bk_results.csv", sep="\t")
	self.df.columns = ['kuta', 'y', 'vr', 'vfr', 'prev']
        # self.df = self.df.drop(self.df[self.df.kuta == "kuta"].index)

        # converting dataframe element types
        self.df["kuta"] = self.df["kuta"].astype('int')
        self.df["y"] = (self.df["y"].astype('float32')).round(decimals=1)
        self.df["vr"] = self.df["vr"].astype('float64')
        self.df["vfr"] = self.df["vfr"].astype('float64')


        self.r = np.unique(self.df.vr.to_numpy(dtype='float64'))

        self.interp_y = []  # to be filled by interpolated functions over r at some fixed y
        for i in np.arange(-4.8, 30.0, 0.1):
            sub_ = self.df.loc[(self.df['kuta'] == 4) & (self.df['y'] == i.round(decimals=1))]
            x_ = np.concatenate(sub_[['vr']].to_numpy(), axis=0)
            y_ = np.concatenate(sub_[['vfr']].to_numpy(), axis=0)

            self.interp_y.append(interpolate.CubicSpline(x_, y_))

        self.interp_r = []  # to be filled by interpolated functions over y at some fixed r
        for r_ in self.r:
            sub_ = self.df.loc[(self.df['kuta'] == 4) & (self.df['vr'] == r_)].sort_values(by='y')
            # sub_ = sub_.sort_values(by='y')
            x_ = np.concatenate(sub_[['y']].to_numpy(), axis=0)
            y_ = np.concatenate(sub_[['vfr']].to_numpy(), axis=0)
            self.interp_r.append(interpolate.CubicSpline(x_, y_))

    def bk_f(self, y):  # fundamental representation

        #        t1 = time.time()
        rap = np.arange(-4.8, 30.0, 0.1)
        # search file to see if interpolation already exists !!!
	print("bk_f y: " + str(y))
        if y in rap:
            i = np.where(rap == y)[0][0]
            return self.interp_y[i]
        else:
            a = np.zeros(len(self.r))
            for i in range(len(a)):
                f = self.interp_r[i]
                a[i] = f(y)
            return interpolate.CubicSpline(self.r, a)

        # write interpolated object to file !!!

    def bk_a(self, y):
	print("bk_a : " + str(y))
        f = self.bk_f(y)
        x = self.r

        g = 2 * f(x) - np.power(f(x), 2)
        return interpolate.CubicSpline(self.r, g)

    def udg_f(self, x, k):
        f = self.bk_f(np.log(self.x0 / x))
        integrand = lambda r: (1 - f(r)) * self.bessel(k * r, 0)
        return 2 * np.pi * intg.quad(integrand, self.xr1, self.xr2, epsabs=1.e-4)[0]

    def udg_a(self, x, k):
        f = self.bk_a(np.log(self.x0 / x))
        integrand = lambda r_: (1 - f(r_)) * self.bessel(k * r_, 0)
        return 2 * np.pi * intg.quad(integrand, self.xr1, self.xr2, epsabs=1.e-4)[0]

    def bessel(self, x, alpha):
        f = lambda t: np.cos(alpha * t - x * np.sin(t))
        return (1 / np.pi) * intg.quad(f, 0, np.pi, epsabs=1.e-3)[0]


# end of class

if __name__ == "__main__":
    n = N()
    x = 0.002463981259066879
    k = 71.1797503982
    x0 = 0.02

    n.bk_f(-1.0)
