import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import csv
import subprocess

import lhapdf as pdf

# flavors: [-5 -4 -3 -2 -1 1 2 3 4 5 21] --> [bb, cb, sb, ub, db, d, u, s, c, b, g]
bb, cb, sb, ub, db, d, u, s, c, b, g = -5, -4 , -3 , -2, -1, 1, 2, 3, 4, 5, 21
ff = pdf.mkPDF("DSS07PI",0)
print(ff.flavors())

qsq2 = 10.  # GeV^2

x = np.linspace(0, 1., num = 1000)
y = [ff.xfxQ2(c, x_, qsq2) + ff.xfxQ2(cb, x_, qsq2) for x_ in x]

plt.xlabel('x')
plt.xscale('log')
# plt.yscale('log')
plt.plot(x,y)
plt.show()





