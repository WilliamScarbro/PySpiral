from mats import *

def reverse(n,r):
    p=getMod(n,NW=True)
    print(p)
    m=n//r
    pr=PolyRing(r,p)
    pr2=PolyRing(m,p)
    t1=Tensor(NTT(pr),Ident(m,p))
    t2=Tensor(Ident(m,p),NTT(pr))
    l1 = LPerm(n,p,r)
    l2 = LPerm(n,p,m)
    print(t2.getMat()) #(NTT_r X I_m)
    print(MM(l2,MM(t1,l1)).getMat()) #L_r (I_m X NTT_r) L_m

def transform(n,r):
    p=getMod(n,NW=True)
    print(p)
    m=n//r
    pr=PolyRing(r,p)
    pr2=PolyRing(m,p)
    t1=Tensor(Ident(r,p),NTT(pr))
    t2=Tensor(NTT(pr2),Ident(r,p))
    l1 = LPerm(n,p,r)
    l2 = LPerm(n,p,m)
    print(t1.getMat()) #(NTT_r X I_m)
    print(MM(l1,MM(t2,l2)).getMat()) #L_r (I_m X NTT_r) L_m


def twiddle(n,r):
    p=getMod(n,NW=False)
    print(p)
    m=n//r
    pr=PolyRing(n,p)
    m1=LPerm(n,p,r)
    m2=Tw_CT(pr,r)
    m3=LPerm(n,p,m)
    print(m1.getMat())
    print(m2.getMat())
    print(m3.getMat())
    print(MM(m1,MM(m2,m3)).getMat())
    print(Tw_CT(pr,m))

def ntt1(S,m):
    p=getMod(S,NW=False)
    n=S//m
    prN=PolyRing(n,p)
    prM=PolyRing(m,p)
    NTTn = NTT(prN)
    NTTm = NTT(prM)
    In = Ident(n,p)
    Im = Ident(m,p)
    Tn = Tw_CT(n,p)
   # dit = MM(MM(Tensor(NTTm,Ident(
#transform(6,3)
#reverse(6,3)
#twiddle(8,4)

pr = PolyRing(4,17)
ntt = NTT(pr)
print(ntt)
gs = GS(ntt,2)
print(gs)

ct = CT(ntt,2)
print(ct)
