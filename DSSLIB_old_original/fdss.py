import dssfh.so as f

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
def fDSS(ih,ic,io,x,q2):
    # return array type [u,ub,d,db,s,sb,c,b,gl]

    [u, ub, d, db, s, sb, c, cb, b, bb, gl] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    f.fDSS(ih,ic,io,x,q2,gl,u,ub,d,db,s,sb,c,cb,b,bb)

    return [u, ub, d, db, s, sb, c, cb, b, bb, gl]