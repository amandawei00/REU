import lhapdf
import numpy as np
import matplotlib.pyplot as plt

# class PDFread():
#    def __init__(self):
#        p = lhapdf.mkPDF("CT10",0)
#        p = lhapdf.mkPDF("C10nlo/0")
#        p.flavors()

#    def pdf(self,flavor,x,q2):
#        return p.xfxQ2(flavor,x,q2)

    
p = lhapdf.mkPDF("CT10",0)
p = lhapdf.mkPDF("CT10/0")
p.flavors()

q2 = 100.0
x = np.arange(0,1.0,0.01)
y = [p.xfxQ2(21,x[i],q2) for i in range(len(x))]
for i in y:
    print(y)

plt.plot(x,y)
plt.yscale('log')
plt.show()





