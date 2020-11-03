import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def equal_arr(a1, a2):
    a1.sort()
    a2.sort()

    for i in range(0, len(a1)):
        if a1[i] != a2[i]:
            # print(str(a1[i]) + ", " + str(a2[i]))
            return False

    return True


df = pd.read_csv("full_bk_results.csv", sep="\t")
df.columns = ['kuta', 'y', 'vr', 'vfr', 'prev']
# df = df.drop(df[df.kuta == "kuta"].index)

# converting dataframe element types
df["kuta"] = df["kuta"].astype('int')
df["y"] = (df["y"].astype('float32')).round(decimals=1)
df["vr"] = df["vr"].astype('float64')
df["vfr"] = df["vfr"].astype('float64')

rap = np.around(np.arange(-4.9, 30.0, 0.1), decimals=2)
vr_original = np.array(df.loc[(df['kuta'] == 4) & (df['y'] == 0.)][['vr']])

for i in range(len(rap)):
    sub_ = df.loc[(df['kuta'] == 4) & (df['y'] == rap[i])]

    vr = np.array(sub_[['vr']])
    vfr = np.array(sub_[['vfr']])
    """print(rap[i])
    print(len(vr))
    print(len(vfr))
    print(len(vr_original))"""

    # if not equal_arr(vr, vr_original):
    #     print("i: " + str(rap[i]))

    plt.plot(vr_original, vfr)

plt.xscale('log')
plt.show()