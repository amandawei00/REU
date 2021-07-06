# sys.path.append("Ogata/python")
# from FBT import FBT
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
import scipy.integrate as intg

# read results.csv file from BK solution to pandas dataframe
df1 = pd.read_csv("results0--5.csv", sep="\s+") # with negative rapidity, needs to be inverted
df2 = pd.read_csv("results0-30.csv", sep="\s+")

df1.columns = ['kuta', 'y', 'vr', 'vfr', 'prev']
df2.columns = ['kuta', 'y', 'vr', 'vfr', 'prev']

df1 = df1.drop(df1[df1.kuta == "kuta"].index)
df2 = df2.drop(df2[df2.kuta == "kuta"].index)

df1["kuta"] = df1["kuta"].astype('int')
df1["y"] = (df1["y"].astype('float32')).round(decimals=1)
df1["vr"] = df1["vr"].astype('float64')
df1["vfr"] = df1["vfr"].astype('float64')

df1 = df1.drop(df1[df1.y == 0.0].index)
df1_ = pd.DataFrame()

for i in np.arange(-4.9, -0.1, 0.1):
    temp = df1[df1["y"] == i]
    df1_ = pd.concat([df1_, temp])
# y_arr = np.arange(-4.9, 0.1, 0.1)
# df1_arr = []

# for i in y_arr:
#     df1_arr.append(df1.loc[df1['y'] == i])

# df = pd.concat(df1_arr)

# df_f = pd.read_csv("results-f.csv", sep="\s+")
# df_f.columns = ['kuta','y','vr','vfr','prev']
# df_f = df_f.drop(df_f[df_f.kuta == "kuta"].index)

df2["kuta"] = df2["kuta"].astype('int')
df2["y"] = (df2["y"].astype('float32')).round(decimals=1)
df2["vr"] = df2["vr"].astype('float64')
df2["vfr"] = df2["vfr"].astype('float64')

dfs = [df1_, df2]

final_df = pd.concat(dfs)
print(final_df)
# write new concatenated dataframe to csv file
final_df.to_csv('full_bk_results-x0-002.csv', sep='\t', header=False, index=False)

