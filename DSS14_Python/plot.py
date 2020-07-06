import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
from matplotlib import rc
rc('text',usetex=True)


df_DSS14_10 = pd.read_csv('output.csv',delim_whitespace = True)

fig = plt.figure()
ax = plt.subplot(121)

ax.plot(df_DSS14_10['z'],df_DSS14_10['z']*(df_DSS14_10['U']+df_DSS14_10['UB']),ls = '-.' , label = r'$\rm u+\bar{u}$')
ax.plot(df_DSS14_10['z'],df_DSS14_10['z']*(df_DSS14_10['UB']),ls = '--' , label = r'$\rm \bar{u}$')
plt.ylabel(r'$\rm z D_{q}^{\pi^+}(z, Q^2 =10)$')
plt.xlabel(r'$\rm z $')
ax.legend()
ax.set_ylim(0,2.1)

ax = plt.subplot(122)
ax.plot(df_DSS14_10['z'],df_DSS14_10['z']*(df_DSS14_10['D']+df_DSS14_10['DB']),ls = '-.' , label = r'$\rm d+\bar{d}$')
ax.plot(df_DSS14_10['z'],df_DSS14_10['z']*(df_DSS14_10['S']),ls = '--' , label = r'$\rm s$')
plt.xlabel(r'$\rm z $')
ax.legend()
ax.set_ylim(0,2.1)
ax.set_yticklabels([])

plt.tight_layout()
fig.savefig('DSS14_10.pdf')
