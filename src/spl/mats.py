from IPython.core.debugger import set_trace
import numpy as np
from poly_ring import *
from primitives import *
import random

# This program mimics the MatSPL representation in Spiral
# General Design:
#     Everything is a Mat
#     Mats are evaluated using getMat, which is done as late as possible to to maximize flexibility
#     There are 4 base Mats (NTT,TW,L,Ident) which correspond to the Spiral definitions
#     There are 2 binary Mat operations (MM (Matrix Multiply), Tensor) which are themselves Mats

class Mat:
    def __init__(self,n,p):
        self.n=n
        self.p=p
        self.isBase=True
    #returns nxn matrix
    def getMat(self):
        pass
    def __str__(self):
        return "{}\n".format(str(self.getMat()))

#provides interface for codegen 
class FullMat(Mat):
    def __init__(self,n,p):
        Mat.__init__(self,n,p)
    def getOutputs(self,inputs,outputs,localVals):
        for i in range(self.n):
            operands=[]
            for j in range(self.n):
                operands.append(Expr("mult",[self.el(i,j),inputs[j]]))
            outputs[i].definition=Expr("add",operands)
        return inputs,outputs,localVals
    def el(i,j):
        pass

class NTT(FullMat):
    def __init__(self,pr,Inv=False,Tr=False):
        FullMat.__init__(self,pr.n,pr.p)
        self.Inv=Inv
        self.Tr=Tr
        self.NW=pr.Pos
        self.pr=pr
        self.p=pr.p
        self.n=pr.n
        self.pr=pr
        self.psi,self.w=pr.getRoots(Inv)
    def getMat(self):
        return np.array([[self.el(i,j) for i in range(self.n)] for j in range(self.n)])
    def el(self,i,j):
        ret=pow(self.w,i*j,self.p)
        if self.Inv:
            i,j=j,i
            ret *= self.pr.inverse(self.n)
        if self.NW:
            ret*=pow(self.psi,j,self.p)
        return ret%self.p
    def symbolic(self,off):
        return "  "*off+f"NTT: n={self.n} p={self.p}"

class DiagMat(Mat):
    def __init__(self,n,p):
        Mat.__init__(self,n,p)
    def getOutputs(self,inputs,outputs,localVals):
        for i in range(self.n):
            outputs[i].definition=Expr("mult",[self.el(i),inputs[i]])
        return inputs,outputs,localVals
        
#used as superclass for all diagonals (as codegen interface)
class Tw(DiagMat):
    def __init__(self,pr,m,Inv=False):
        DiagMat.__init__(self,pr.n,pr.p)
        self.m=m
        self.psi,self.w=pr.getRoots(Inv)
    def func(self,i,j):
        pass
    def el(self,k):
        return self.func(k//self.m,k%self.m)
    def diagonal(self):
        vec=[]
        for i in range(self.n//self.m):
            for j in range(self.m):
                vec.append(self.func(i,j))
        return vec 
    def getMat(self):
        vec=self.diagonal()
        ret = np.array([[vec[i] if i==j else 0 for i in range(self.n)] for j in range(self.n)])
        return ret
            
    
class Tw_CT(Tw):
    def __init__(self,pr,m,Inv=False):
        Tw.__init__(self,pr,m,Inv)
        self.NW=pr.Pos
        self.n=pr.n
    def func(self,i,j):
         if self.NW:
            return pow(self.psi,(2*i*j+(1)*i)%(2*self.n),self.p)#return (pow(self.w,i*j,self.p)*pow(self.psi,((1-self.m)*i)%(2*n),self.p))%self.p
         else:
            return pow(self.w,i*j,self.p)
        
class Tw_GS(Tw):
    def __init__(self,pr,m,Inv=False):
        Tw.__init__(self,pr,m,Inv)
        self.NW=pr.Pos
    def func(self,i,j):
        if self.NW:
            return pow(self.psi,i+2*i*j,self.p)
        else:
            return pow(self.w,i*j,self.p)
        
class Tensor(Mat):
    def __init__(self,m1,m2):
        Mat.__init__(self,m1.n*m2.n,m1.p)
        assert m1.p==m2.p
        self.m1=m1
        self.m2=m2
        self.isBase=False
    def getMat(self):
        m1m=self.m1.getMat()
        m2m=self.m2.getMat()
        ret = [[0 for i in range(self.n)] for j in range(self.n)]
        for i in range(self.n):
            for j in range(self.n):
                ret[i][j]=m1m[i//m2m.shape[0]][j//m2m.shape[0]]*m2m[i%m2m.shape[0]][j%m2m.shape[0]]%self.p
        return np.array(ret)
    def getOutputs(self,inputs,outputs,localVals): 
        if isinstance(self.m1,DiagMat):
            for i in range(self.m1.n):
                inputs[self.m2.n*i:self.m2.n*(i+1)],outputs[self.m2.n*i:self.m2.n*(i+1)],localVals = self.m2.getOutputs(inputs[self.m2.n*i:self.m2.n*(i+1)],outputs[self.m2.n*i:self.m2.n*(i+1)],localVals)
                for j in range(self.m2.n):
                    #l = Val("local",definition=outputs[self.m2.n*i+j].definition,local=True)
                    #localVals.append(l)
                    outputs[self.m2.n*i+j].definition=Expr("mult",(self.m1.el(i),outputs[self.m2.n*i+j].definition))
            return inputs,outputs,localVals
        else:
            if not isinstance(self.m2,DiagMat):
                print("Error: One operand of tensor must be diagonal!")
            reverse = MM(MM(LPerm(self.n,self.p,self.m1.n),Tensor(self.m2,self.m1)),LPerm(self.n,self.p,self.m2.n))
            return reverse.getOutputs(inputs,outputs,localVals)

class MM(Mat):
    def __init__(self,m1,m2):
        Mat.__init__(self,m1.n,m1.p)
        assert m1.n==m2.n and m1.p==m2.p
        self.m1=m1
        self.m2=m2
    def getMat(self):
        ret = self.m1.getMat()@self.m2.getMat()
        return np.array([[ret[j][i]%self.p for i in range(self.n)] for j in range(self.n)])
    def getOutputs(self,inputs,outputs,localVals):
        if isinstance(self.m2,PermMat):
            inputs = [inputs[self.m2.perm(i)] for i in range(len(inputs))] 
            return self.m1.getOutputs(inputs,outputs,localVals)
        if isinstance(self.m1,PermMat):
            outputs = [outputs[self.m1.iperm(i)] for i in range(len(outputs))]
            return self.m2.getOutputs(inputs,outputs,localVals)
        inputs,outputs,localVals = self.m2.getOutputs(inputs,outputs,localVals)
        interVals=[Val("local",o.definition,local=True) for o in outputs]
        localVals+=interVals
        return self.m1.getOutputs(interVals,outputs,localVals)

class Ident(DiagMat):
    def __init__(self,n,p):
        DiagMat.__init__(self,n,p)
    def el(self,i):
        return 1
    def getMat(self):
        return np.array([[1 if i==j else 0 for i in range(self.n)] for j in range(self.n)])

class PermMat(Mat):
    def __init__(self,n,p):
        Mat.__init__(self,n,p)
    def getOutputs(self,inputs,outputs,localVals):
        for i in range(self.n):
            outputs[i].definition=inputs[self.perm(i)]
        return inputs,outputs,localVals

class LPerm(PermMat):
    def __init__(self,n,p,stride):
        PermMat.__init__(self,n,p)
        self.str=stride
        self.m=n//stride
    def perm(self,i):
        return i//self.m+self.str*(i%self.m)
    def iperm(self,i):
        return i//self.str+self.m*(i%self.str)
    def getMat(self):
        return np.array([[1 if self.perm(j)==i else 0 for i in range(self.n)] for j in range(self.n)])
 
class Reduction(Mat):
    def __init__(self,mats,par):
        super().__init__(par.n,par.p)
        self.isBase=False
        self.mats=mats
        self.par=par
    def multiply(self,x):
        for i in range(len(self.mats)-1,-1,-1):
            x=self.mats[i].getMat()@x
            for j in range(len(x)):
                x[j]%=self.p
        return x   
    def result(self):
        if len(self.mats)==0:
            return "no mats"
        ret=self.mats[0]
        for i in range(1,len(self.mats)):
            ret=MM(ret,self.mats[i]) #removed ret".copy()"
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
    def getOutputs(self,inputs,outputs,localVals):
        return self.result().getOutputs(inputs,outputs,localVals)


class CT(Reduction):
    def __init__(self,ntt,radix,omit_perm=False):
        Reduction.__init__(self,[],ntt)
        #assert ntt.Tr==False
        m=ntt.n//radix
        self.radix=radix
        self.c1 = NTT(ntt.pr.toNeg().reduce(m),Inv=ntt.Inv)
        self.c2 = NTT(ntt.pr.reduce(radix),Inv=ntt.Inv)
        m1=Tensor(self.c1,Ident(m,ntt.p))
        m2=Tw_CT(ntt.pr,m,Inv=ntt.Inv)
        m3=Tensor(Ident(radix,ntt.p),self.c2)
        m4=LPerm(ntt.n,ntt.p,radix)
        self.mats=[m1,m2,m3] if omit_perm else [m1,m2,m3,m4]
    def as_func(self):
        return f"CT: n={self.n} radix={self.radix}"
    
class GS(Reduction):
    def __init__(self,ntt,radix,omit_perm=False):
        Reduction.__init__(self,[],ntt)
        #assert ntt.Tr==True
        m=ntt.n//radix
        self.radix=radix
        self.c1 = NTT(ntt.pr.reduce(m),Inv=ntt.Inv,Tr=True)
        self.c2 = NTT(ntt.pr.toNeg().reduce(radix),Inv=ntt.Inv,Tr=True)
        m1=LPerm(ntt.n,ntt.p,radix)
        m2=Tensor(Ident(m,ntt.p),self.c1)
        m3=Tw_CT(ntt.pr,radix,Inv=ntt.Inv)
        m4=Tensor(self.c2,Ident(radix,ntt.p))
        self.mats=[m2,m3,m4] if omit_perm else [m1,m2,m3,m4]
    def as_func(self):
        return f"GS: n={self.n} radix={self.radix}"

class Spiral_DIF(Reduction):
    def __init__(self,ntt,radix,omit_perm=False):
        Reduction.__init__(self,[],ntt)
        #assert ntt.Tr==True
        m=ntt.n/radix
        self.c1 = NTT(ntt.pr.toNeg().reduce(m),Inv=ntt.Inv)
        self.c2 = NTT(ntt.pr.reduce(radix),Inv=ntt.Inv)
        m1=Tensor(self.c1,Ident(m,ntt.p))
        m2=Tw_CT(ntt.pr,radix,Inv=ntt.Inv)
        m3=Tensor(Ident(radix,ntt.p),self.c2)
        m4=LPerm(ntt.n,ntt.p,radix)
        self.mats=[m1,m2,m3] if omit_perm else [m1,m2,m3,m4]
    def as_func(self):
        return f"CT: n={self.n} radix={self.radix}"
    
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
