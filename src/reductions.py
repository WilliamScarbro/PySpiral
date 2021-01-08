#functions to define and test transforms
import random
from mats import *
from poly_ring import *
from mats import *

class Reduction(Mat):
    def __init__(self,mats,par):
        super().__init__(par.n,par.p)
        self.isBase=False
        self.mats=mats
        self.par=par
    def multiply(self,x):
        for i in range(len(self.mats)-1,-1,-1):
            print(x)
            x=self.mats[i].getMat()@x
            for j in range(len(x)):
                x[j]%=self.p
        return x   
    def result(self):
        if len(self.mats)==0:
            return "no mats"
        ret=self.mats[0]
        for i in range(1,len(self.mats)):
            ret=MM(ret.copy(),self.mats[i])
        return ret
    def __str__(self):
        ret=""
        for m in self.mats:
            ret+=str(m)
        return ret
    def getMat(self):
        t_res=self.mats[0].getMat()
        for i in range(1,len(self.mats)):
            t_res=t_res@self.mats[i].getMat()
        t_res%=self.p
        return t_res
    def test_transform(self):
        t_res=self.getMat()
        print("reduction:\n%s"%t_res)
        print("parent:\n%s"%self.par.getMat())
        return t_res==self.par.getMat()
    def symbolic(self,off):
        ret = "  "*off+self.as_func()
        ret += "\n"+self.c1.symbolic(off+1)
        ret += "\n"+self.c2.symbolic(off+1)
        return ret
class CT(Reduction):
    def __init__(self,ntt,radix,omit_perm=False):
        Reduction.__init__(self,[],ntt)
        assert ntt.Tr==False
        m=ntt.n//radix
        self.radix=radix
        self.c1 = NTT(ntt.pr.toNeg().reduce(m),Inv=ntt.Inv)
        self.c2 = NTT(ntt.pr.reduce(radix),Inv=ntt.Inv)
        m1=Tensor(self.c1,I(m,ntt.p))
        m2=Tw_CT(ntt.pr,m,Inv=ntt.Inv)
        m3=Tensor(I(radix,ntt.p),self.c2)
        m4=L(ntt.n,ntt.p,radix)
        self.mats=[m1,m2,m3] if omit_perm else [m1,m2,m3,m4]
    def as_func(self):
        return f"CT: n={self.n} radix={self.radix}"
    
class GS(Reduction):
    def __init__(self,ntt,radix,omit_perm=False):
        Reduction.__init__(self,[],ntt)
        assert ntt.Tr==True
        m=ntt.n//radix
        self.radix=radix
        self.c1 = NTT(ntt.pr.reduce(m),Inv=ntt.Inv,Tr=True)
        self.c2 = NTT(ntt.pr.toNeg().reduce(radix),Inv=ntt.Inv,Tr=True)
        m1=L(ntt.n,ntt.p,radix)
        m2=Tensor(I(m,ntt.p),self.c1)
        m3=Tw_CT(ntt.pr,radix,Inv=ntt.Inv)
        m4=Tensor(self.c2,I(radix,ntt.p))
        self.mats=[m2,m3,m4] if omit_perm else [m1,m2,m3,m4]
    def as_func(self):
        return f"GS: n={self.n} radix={self.radix}"
    
def convolution(pr,x,y):
    ret=[]
    for i in range(pr.n):
        a=sum(x[j]*y[i-j] for j in range(i))
        b=sum(x[j]*y[pr.n+i-j] for j in range(i+1,pr.n))
        c = (a-b)%pr.p if pr.Pos else (a+b)%pr.p
        ret.append(c)
    return np.array(ret)

#reads radices as binary heap, 
#  radix 1 represents no change, 
#  numbers appearing below a 1 will not be read
#returns reduced 
def build_reduction(cur,red,radices,i,omit_perm=False):
    if i>=len(radices) or radices[i]==1:
        return cur
    assert cur.n>2, f"Applying reduction to NTT_{cur.n}<=2 not allowed"
    reduced = red(cur,radices[i],omit_perm=omit_perm)
    reduced.c1 = build_reduction(reduced.c1,red,radices,2*i+1,omit_perm=omit_perm)
    reduced.c2 = build_reduction(reduced.c2,red,radices,2*i+2,omit_perm=omit_perm)
    return reduced

def build_reduction_tt(cur,red,radices,omit_perm=False):
    if isinstance(radices,int):
        return red(cur,radices,omit_perm=omit_perm) if radices!=1 else cur
    assert len(radices)==2
    assert len(radices[1])==2
    radix = radices[0]
    lradices = radices[1][0]
    rradices = radices[1][1]
    assert cur.n>2, f"Applying reduction to NTT_{cur.n}<=2 not allowed"
    reduced = red(cur,radix,omit_perm=omit_perm)
    #the children are reversed between CT vs GS, but this is not represented here...
    reduced.c1 = build_reduction_tt(reduced.c1,red,lradices,omit_perm=omit_perm)
    reduced.c2 = build_reduction_tt(reduced.c2,red,rradices,omit_perm=omit_perm)
    return reduced

#converts binary tree (list) representation to tuple tree representation
def bt2tt(bt,i):
    if 2*i+1>=len(bt) or bt[i]==1:
        return bt[i]
    left = bt2tt(bt,2*i+1)
    right = bt2tt(bt,2*i+2) if 2*i+2<len(bt) else 1
    return (bt[i],(left,right))

def random_tt(n):
    factors=[]
    for i in range(2,int(n**1/2)+1):
        if n%i==0:
            factors.append(i)
    if not factors:
        return n
    radix=random.choice(factors)
    return (radix,(random_tt(radix),random_tt(n//radix))) if radix!=2 else 2

#changes radices to produce the inverse permutation
#radices must be tuple tree
def convert_radices(n,radices):
    if isinstance(radices,int):
        return n//radices
    return (n//radices[0],(radices[1][1],radices[1][0]))
    
def Mv(pr,m,v):
    return np.mod(m.getMat()@v,pr.p)
