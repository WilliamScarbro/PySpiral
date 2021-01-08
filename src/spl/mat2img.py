from mats import *
from poly_ring import *
from PIL import Image
import numpy as np


def el2color(el):
    if el==0:
        return (255,255,255)
    val = hash((el+0.1)**4)
    #print(val)
    return (val%256,val/256%256,val/256/256%256)

def mat2pix(mat,shape):
    pixels=mat.getMat()
    print(pixels)
    rep=shape//pixels.shape[0]
    pixels=np.array([[el2color(pixels[i//rep][j//rep]) for i in range(shape)] for j in range(shape)],dtype=np.uint8)
    return pixels
 
def mat2img(mat,name,shape):
    # Use PIL to create an image from the new array of pixels
    pixels=mat2pix(mat,shape)
    new_image = Image.fromarray(pixels)
    new_image.save(name+".png")

def red2img(red,name,shape):
    i=0
    for m in red.mats:
        pixels=mat2pix(m,shape)
        new_image = Image.fromarray(pixels)
        new_image.save(name+f"{i}.png")
        i+=1
        
n=16
p=getMod(16,NW=True)
size = 400
r=4

pr=PolyRing(n,p,Pos=True)
ntt=NTT(pr)
mat2img(ntt,"ntt",size)

i = Ident(n,p)
mat2img(i,"ident",size)


l = LPerm(n,p,r)
mat2img(l,"lPerm",size)

ct = CT(ntt,r)
red2img(ct,"CT",400)
