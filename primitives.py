class Expr:
    def __init__(self,operation,operands):
        self.operation=operation
        self.operands=operands
    def __eq__(self,obj):
        if not isinstance(obj,Expr):
            return False
        if self.operation.commutative:
            pass
            #combare sets of operands
        else:
            pass
            #compare order of operands
    def __str__(self):
        ret=f"{self.operation}("
        for op in self.operands:
            ret+=str(op)
        return ret+")" 

#String: label is name of variable
#Expr: definition defines how it is calculated (blank for inputs)
class Val(Expr):
    def __init__(self,label,definition=None,local=False):
        self.label=label
        self.definition=definition
    def __str__(self):
        return str(self.label)
    def toSymbolic(self,interVals):
        return [Expr("assign",[self,visit(self.definition)])]

#expr can be Expr,Val,str,int
def visit(expr):
    if isinstance(expr,Val):
        return expr
    if isinstance(expr,Expr):
        if expr.operation=="add":
            newOps=[]
            for op in expr.operands:
                if op!=0:
                    newOps.append(visit(op))
            expr.operands=newOps
            return expr
        if expr.operation=="mult":
            newOps=[]
            for op in expr.operands:
                if op!=1:
                    newOps.append(visit(op))
            expr.operands=newOps
            return expr
    return expr
