import numpy as np
import scipy.interpolate as interpolate
import pandas as pd
import csv

class N():
    def __init__(self):
        n = 400

# read results.csv file from BK solution to pandas dataframe
        df = pd.read_csv("results.csv", sep="\s+")
        df = df.drop(df[df.kuta == "kuta"].index)

# converting dataframe element types
        df["kuta"] = df["kuta"].astype('int')
        df["y"] = (df["y"].astype('float32')).round(decimals=1)
        df["vr"] = df["vr"].astype('float64')
        df["vfr"] = df["vfr"].astype('float64')

# isolating relevant data
        df = df.drop(df[(df.y != y) | (df.kuta != 4)].index)

# assigning vr and vfr values to xlist, ylist, for interpolation
        r = df.vr.to_numpy(dtype='float64')
        y = df.y.to_numpy(dtype='float32')
        n = df.vfr.to_numpy(dtype='float64')

        self.f = interpolate.interp2d(r, y, n, kind=3) #WRITE METHOD TO FIND FOURIER TRANSFORM OF N, AND F(A)

    # fundamental representation
    def bk_f(self,x2,y): 
        # y: rapidity
        return self.f(x2,y)

    def bk_a(self,x2,y):
        n_a = 2*self.bk_f(x2,y) - np.power(self.bk_f(x2,y),2)
        return n_a


