import dssfh as f
import matplotlib.pyplot as plt
import numpy as np
"""  
    input: 
        ih -> hadron type (1: pion, 2: kaon, 3: proton, 4: charged hadrons)
        ic -> hadron charge (0: 0, 1: +, -1: -)
        io -> order (0: LO, 1: NLO)
        x -> 0.05 < x < 1.0
        q2 -> scale in GeV^2, between 1.0 and 1.e5
    output:
        u
        ub
        d 
        db
        s 
        sb
        c
        b
        gl
"""
print(f.__doc__)
def ff(ih,ic,io,x,q2):
    # return array type [u,ub,d,db,s,sb,c,b,gl]

#    [u, ub, d, db, s, sb, c, cb, b, bb, gl] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
     gl,u,ub,d,db,s,sb,c,cb,b,bb = f.fdss(ih,ic,io,x,q2)

     return u

x_=np.arange(0.5,1,10)
ih=2
ic=1
io=0
q2=10

print(ff(ih,ic,io,0.9,q2))
#y_=[ff(ih,ic,io,x,q2)[gl] for x in x_]

#plt.plot(x,y)
#plt.show()
