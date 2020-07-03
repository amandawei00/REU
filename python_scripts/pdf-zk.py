
import numpy as np
import matplotlib.pyplot as plt

# for details, see: https://arxiv.org/pdf/1412.7420.pdf

# PDFs: parton distribution function
# both call the central number (0) of CT14nlo
p = lhapdf.mkPDF("CT14nlo", 0)
p = lhapdf.mkPDF("CT14nlo/0")

# let us look at the flavors
# they follow convention of particle data group
# d = 1, u = 2, s = 3, c = 4, b = 5, g = 21
p.flavors()

x = 0.5
Q2 = 100.0
d = p.xfxQ2(1, x, Q2)
u = p.xfxQ2(2, x, Q2)
s = p.xfxQ2(3, x, Q2)
c = p.xfxQ2(4, x, Q2)
b = p.xfxQ2(5, x, Q2)
g = p.xfxQ2(21, x, Q2)

print('%0.2e, %0.2e, %0.2e, %0.2e, %0.2e, %0.2e' %(d,u,s,c,b,g))

ff = lhapdf.mkPDF("JAM19FF_pion_nlo", 0)
df = ff.xfxQ2(1, x, Q2)
uf = ff.xfxQ2(2, x, Q2)
sf = ff.xfxQ2(3, x, Q2)
cf = ff.xfxQ2(4, x, Q2)
bf = ff.xfxQ2(5, x, Q2)
gf = ff.xfxQ2(21, x, Q2)
print('%0.2e, %0.2e, %0.2e, %0.2e, %0.2e, %0.2e' %(df,uf,sf,cf,bf,gf))

x = np.empty(100)
y = np.empty(100)
z = np.empty(100)
for i in range(100):
    x[i] = 0.01 * i + 0.01
    y[i] = p.xfxQ2(21, x[i], Q2)
    z[i] = ff.xfxQ2(21, x[i], Q2)

plt.plot(x, y, 'ro', label='g-PDF', markersize=3)
plt.plot(x, z, 'b^', label='g-FF', markersize=3)
plt.yscale('log')
# plt.xscale('log')
plt.legend()
plt.show()










