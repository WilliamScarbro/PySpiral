import sys
sys.path.insert(1,"../spl")
from mats import *
from poly_ring import *


def test_ct():
    pr4 = PolyRing(4)
    pr2 = PolyRing(2,17)
    ntt4 = NTT(pr4)
    ntt2 = NTT(pr2)
    ct = CT(ntt4,2)
    print(ntt4.getMat())
    print(ct)
test_ct()
