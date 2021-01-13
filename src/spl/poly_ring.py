

#Contains PolyRing and other ff related functions

class PolyRing:
    def __init__(self,n,p=None,Pos=False):
        self.n=n
        self.p=p if p!=None else getMod(self.n,NW=not Pos)
        self.Pos=Pos
        assert (self.p-1)%n==0
        k=(self.p-1)//n
        r=primitive_root(self.p)
        self.psi=pow(r,k//2,self.p)
        self.w=pow(r,k,self.p)
        self.psi_i = pow(self.psi,2*n-1,self.p)
        self.w_i = pow(self.w,n-1,self.p)
        assert (self.psi*self.psi_i) % self.p == 1
        assert (self.w*self.w_i) % self.p == 1
    def __str__(self):
        return f"PolyRing Z/{self.p}[X]/(X^{self.n} {'+' if self.Pos else '-'} 1)\n "
    def getRoots(self,Inv):
        if Inv:
            return self.psi_i,self.w_i
        else:
            return self.psi,self.w
    def reduce(self,m):
        assert self.n%m==0
        return PolyRing(self.n//m,self.p,Pos=self.Pos)
    def toNeg(self):
        return PolyRing(self.n,self.p,Pos=False)
    def inverse(self,a):
        return pow(a,self.p-2,self.p)
    
#some helper functions
def equal(m1,m2):
    if m1.n!=m2.n or m1.p!=m2.p:
        return False
    m1m=m1.getMat()
    m2m=m2.getMat()
    for i in range(m1.n):
        for j in range(m1.n):
            if m1m[i][j]!=m2m[i][j]:
                return False
    return True

def transform(m,v):
    mm=m.getMat()
    ret = mm*v
    return [i%m.p for i in ret]


def primitive_root(p):
    for i in range(1,p-1):
        eq1=False
        a=1
        for j in range(p-2):
            a=(a*i)%p
            if a==1:
                eq1=True
                break
        if not eq1:
            return i
    return -1

def is_prime(x):
    for i in range(2,int(x**1/2)+1):
        if x%i==0:
            return False
    return True

def getMod(n,NW=False):
    if NW:
        n+=n
    c=n
    while (not is_prime(c+1)):
        c+=n
    return c+1
