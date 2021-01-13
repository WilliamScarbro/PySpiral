import sys
sys.path.insert(1,"../spl")
import codegen
import mats
import poly_ring
import numpy
import subprocess
import os


# This suite tests generated c code against the transforms from mats
# This assumes that all of the transforms from mats are correct

def main_for_test(test_location,test_name, inputs, dtype):
    length = len(inputs)
    content = f'#include <stdio.h>\n#include "{test_location}"\n\n'
    content+="int main(int argc, char** argv){\n"
    content+=f"{dtype} X[{length}] = "+'{'+str(inputs)[1:-1]+'};\n'
    content+=f"{dtype} Y[{length}];\n" 
    content+=f"{test_name}(X,Y);\n"
    content+='printf("[%ld",Y[0]);\nfor(int i=1; i<'+str(length)+'; i++){\nprintf(",%ld",Y[i]);\n}\nprintf("]");}'
    return content

def test_transform(mat,name,inputs=None):
    if inputs==None:
        inputs=[i for i in range(mat.n)]
    correct_out = mat.getMat()@numpy.array(inputs).reshape(-1,1) #will change once reductions are implemented
    codegen.codeWrite(mat,name,f"../../bin/{name}.c")
    main_contents = main_for_test(name+".c",name,inputs,"long")
    main_c = open("../../bin/main.c","w")
    main_c.write(main_contents)
    main_c.close()
        
    os.system(f"gcc ../../bin/main.c -o ../../bin/{name}.o")
    out = subprocess.check_output(f"../../bin/{name}.o",shell=True)
    decoded = out.decode("utf-8")
    out = [int(d) for d in decoded[1:-1].split(",")]
    print(out)
    print(correct_out.reshape(1,-1))
    print(numpy.array(out).reshape(-1,1)==correct_out) 

def test_ntt():
    pr = poly_ring.PolyRing(8)
    ntt = mats.NTT(pr)
    test_transform(ntt,"ntt_8")

def test_tensor():
    pr=poly_ring.PolyRing(4)
    ntt4 = mats.NTT(pr)
    t = mats.Tensor(mats.Ident(2,17),ntt4)
    test_transform(t,"tensor_8")

def test_gs():
    pr=poly_ring.PolyRing(4)
    ntt4 = mats.NTT(pr)
    gs=mats.GS(ntt4,2)
    print(gs)
    test_transform(gs,"GS_4")

test_gs()
#test_tensor()
#test_ntt()    
#main_for_test("file_name","test_name",[0,1,2,3,4,5,6,7],"long")
