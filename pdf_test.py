import numpy as np
import matplotlib.pyplot as plt
import scipy

import lhapdf as pdf

bb, cb, sb, ub, db, d, u, s, c, b, g = -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21
pdf = pdf.mkPDF("CT10",0)

f = g  # flavor
q = [10, 50, 1.e2, 2.e2, 5.e2, 1.e3, 2.e3, 5.e3]
qsq2 = [np.power(q[i],2) for i in range(len(q))]

x = np.linspace(1.e-8, 1.e0, 100000000)
y = []

for i in range(len(qsq2)):
	y.append([pdf.xfxQ2(f, x_, qsq2[i]) for x_ in x])

fig, ax = plt.subplots(1,1)
for i in range(len(qsq2)):
	ax.plot(x, y[i], label="qsq2="+str(qsq2[i]))

ax.set_xlim(1.e-8, 1.e0)
ax.set_xlabel('x')
ax.set_xscale('log')


# plt.xlabel('x')
# plt.xscale('log')
# plt.legend()


plt.show()



