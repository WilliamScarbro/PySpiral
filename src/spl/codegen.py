from mats import *
from poly_ring import *
from primitives import *

#class Transfrom:
#    def __init__(self,inputS,outputS,mat):
        
def codeGen(mat,fname):
    #use X as inputs and Y as result
    inputs = [Val(f"X[{i}]") for i in range(mat.n)]
    outputs = [Val(f"Y[{i}]") for i in range(mat.n)]
    localVals = []


    #step 1: build expression tree from mat (recursive if mat is a reduction)
    #each output (layer by layer) has an associated expression tree     
    inputs,outputs,localVals = exprGen(mat,inputs,outputs,localVals)
    #outputs and localVals all be defined

    #step 2: convert expression tree to symbolic code (still contains operation names)  
    #order of computation: inputs -> localVals (in order) -> outputs
    symbolLines, interVals = symbolGen(inputs,outputs,localVals)
    
    #step 3: convert symbolic code to C (start with template, define local varaibles, convert symbolic operations to valid C)
    return cGen(fname,inputs,outputs,localVals+interVals,symbolLines)

#caller is reponsible for creating outputs and adding to localVals
def exprGen(mat,inputs,outputs,localVals):
   return mat.getOutputs(inputs,outputs,localVals)


def symbolGen(inputs,outputs,localVals):
    symbolLines = []
    interVals=[] #need to be initialized
    for lv in localVals:
        symbolLines+=lv.toSymbolic(interVals)
    for out in outputs:
        symbolLines+=out.toSymbolic(interVals)
    return symbolLines, interVals

def cGen(fName,inputs,outputs,localVals,symbolLines):
    from datetime import date
    vType = "long"
    start_template=f"\* sporbiter generated code \n   created {date.today().strftime('%m/%d/%y')}*\ \n"
    function_wrapper=f"void {fName}({vType}* X, {vType}* Y)"+"{\n"
    if localVals:
        init_locals=f"  {vType} {localVals[0].label}"
        for lv in localVals[1:]:
            init_locals+=", "+lv.label
        init_locals+=";\n"
    else:
        init_locals="  \*no locals *\ \n"
    computation=""
    equivMap={}
    for sl in symbolLines:
        code = expr2C(sl)
        if code !="":
            computation+="  "+code+";\n"
    return start_template+function_wrapper+init_locals+computation+"}\n"

#expr can be Expr or int
#all Vals should have been replaced by assign(label,definition)
def expr2C(expr):
    if isinstance(expr,Val):
        #if expr.label in equivMap:
        #    print("replaced equivalence")
        #    return equivMap[expr.label].label
        return expr.label
    if isinstance(expr,str):
        return expr
    if isinstance(expr,int):
        return str(expr)
    if expr.operation=="assign":
        assert isinstance(expr.operands[0],Val)
        if isinstance(expr.operands[1],Val):
            #print(f"found equiv {expr.operands[0]}={expr.operands[1]}")
            #equivMap[expr.operands[0].label]=expr.operands[1].label
            expr.operands[0].setEqual(expr.operands[1])
            return ""
        #print(f"assignment: {expr.operands[1]}")
        return expr.operands[0].label+"="+expr2C(expr.operands[1])
    if expr.operation=="add":
        addition = expr2C(expr.operands[0])
        for op in expr.operands[1:]:
            addition += "+"+expr2C(op)
        return addition
    if expr.operation=="mult":
        multiplication = expr2C(expr.operands[0])
        for op in expr.operands[1:]:
            multiplication += "*"+expr2C(op)
        return multiplication
    return f"Error: no matches found for expr: {expr}"

if __name__=="__main__":


    pr2=PolyRing(2,17)
    ntt2=NTT(pr2,Inv=False)
    print(codeGen(ntt2,"NTT2"))

    pr4=PolyRing(4,17)
    ntt4=NTT(pr4,Inv=False)
    print(codeGen(ntt4,"NTT4"))

    pr8 = PolyRing(8,17)
    ntt8=NTT(pr8,Inv=False)
    
    tw = Tw_CT(pr4,2)
    print(codeGen(tw,"Tw4_2"))

    i = Ident(pr4.n,pr4.p)
    print(codeGen(i,"I4"))

    l = LPerm(4,17,2)
    print(codeGen(l,"L4_2"))

    t = Tensor(Ident(2,17),ntt4)
    print(codeGen(t,"T_4"))

    t2 = Tensor(ntt4,Ident(2,17))
    print(codeGen(t2,"T_4_rev"))

    compose = MM(LPerm(8,17,2),MM(t2,LPerm(8,17,4)))
    print(codeGen(compose,"composed_tensor"))

    m = MM(l,ntt4)
    print(codeGen(m,"MM_4"))

    m2 = MM(ntt4,l)
    print(codeGen(m2,"MM_4_rev"))

    gs = GS(ntt8,2)
    print(codeGen(gs,"GS_NTT8_2"))

    ct = CT(ntt8,2)
    print(codeGen(ct,"CT_NTT8_2"))
