import numpy as np
import scipy.interpolate as interpolate
import pandas as pd
import csv

def bk_func(y):
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
    x = df.vr.to_numpy(dtype='float64')
    y = df.vfr.to_numpy(dtype='float64')

    f = interpolate.interp1d(x, y, kind=5)



bk_func(1.0)




