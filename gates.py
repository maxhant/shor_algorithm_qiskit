import qiskit as q
import numpy as np

def Carry(c, idx, reverse = False):
    """
    Function: 
    It does bit's addition by considering a previously carried bit.
    It return a carry bit.
    
    Args:
    c: circuit
    idx: list of 4 qubits [carry_i, bit_1, bit_2, carry_(i+1)]
    reverse (opt): bool to run the operation in reverse
    
    Return:
    The circuit which carries intrisically the qubit's state. 
    Carry_(i+1) carries the excess forward
    """
    if not len(idx) == 4:
        raise Exception("Carry gate takes 4 qubits")
    if not reverse:
        c.ccx(idx[1],idx[2],idx[3])
        c.cx(idx[1],idx[2])
        c.ccx(idx[0],idx[2],idx[3])
    else:
        c.ccx(idx[0],idx[2],idx[3])
        c.cx(idx[1],idx[2])
        c.ccx(idx[1],idx[2],idx[3])      
    c.to_gate(label="Carry")
    return c

def Sum(c, idx, reverse = False):
    """
    Function: 
    It set the carry qubits back to 0
    
    Args:
    c: circuit
    idx: list of 3 qubits [carry_i, bit_1, bit_2]
    reverse (opt): bool to run the operation in reverse
    
    Return:
    The circuit which carries intrisically the qubit's state. 
    Unchanged bit_1 and bit_2 but carry_i is set back to 0 if bit_1 and bit_2 
    are ones.
    """
    if not len(idx) == 3:
        raise Exception("Carry gate takes 4 qubits")        
    if not reverse:
        c.cx(idx[1], idx[2])
        c.cx(idx[0], idx[2])
    else:
        c.cx(idx[0], idx[2])
        c.cx(idx[1], idx[2])    
    c.to_gate(label="Sum")
    return c

def Adder(c, reg1, reg2, ancil_adder, reverse = False):
    """
    Function: 
    It adds two numbers together. 
    
    Args:
    c: circuit
    reg1: one of the number added of length n
    reg2: the other number of length n
    ancil_adder: the n carrying qubits required
    reverse (opt): bool to run the operation in reverse
    
    Return:
    The circuit which carries intrisically the qubit's state. 
    reg1 stays unchanged.
    reg2 become reg1+reg2 or reg2-reg1 for reverse = False and True respectively
    ancil_adder becomes completely 0 
                except for its first qubit which indicates overflow
    """
    n = len(reg1)
    try:
        ancil = ancil_adder[0:n+1]
    except: 
        raise Exception(f" {len(ancil_adder)} ancillaes when {n+1} required.")
    if n != len(reg2):
        raise Exception("Currently only supports addition of equal bit length")        
    if not reverse:
        for i in range(n):
            i = (n-1) - i
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]])
        c.cx(reg1[0], reg2[0])
        c = Sum(c, [ancil[1], reg1[0], reg2[0]])
        for i in range(1, n):
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]], reverse=True)
            c = Sum(c, [ancil[i+1], reg1[i], reg2[i]])        
    else:
        for i in range(1, n):
            i = n - i
            c = Sum(c, [ancil[i+1], reg1[i], reg2[i]], reverse=True)            
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]], reverse=False)
        c = Sum(c, [ancil[1], reg1[0], reg2[0]], reverse=True)
        c.cx(reg1[0], reg2[0])
        for i in range(n):
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]], reverse=True) 
    c.to_gate(label="Adder")
    return c

def Adder_mod(c, reg1, reg2, regN, regN_ctrl, ancil_mod, reverse=False):
    """
    Function: 
    It takes the modulo of two numbers' sum. 
    
    Args:
    c: circuit
    reg1: one of the number added of length n (reg1 > 0)
    reg2: the other number of length n, note that reg2 < N
    regN: register used to stored N used for the modulo
    regN_ctrl: register used to copy N
    ancil_mod: the n carrying qubits required and one extra controlled qubit 
    reverse (opt): bool to run the operation in reverse
    
    Return:
    The circuit which carries intrisically the qubit's state. 
    reg1 stays unchanged.
    reg2 becomes (reg1+reg2)%N for reverse = False 
    reg2 becomes reg2-reg1 if reg2-reg1>N else N+reg2-reg1 for reverse = True
    ancil_mod in fully set to 0
    reg_N remains unchanged
    reg_N_ctrl is set back to 0
    """
    n = len(reg1)
    ancil = ancil_mod[0:n+1]
    if n != len(reg2):
        raise Exception("Currently only supports addition of equal bit length")
    try:
        ctrl = ancil_mod[n+1]
    except:
        raise Exception(f"Not enough ancilla with {len(ancil_mod)} given for {n+2} required")
    
    if not reverse:
        c = Adder(c, reg1, reg2, ancil, reverse = False)
        c = Adder(c, regN, reg2, ancil, reverse = True)
        c.cx(ancil[0], ctrl)
        for i in range(n):
            c.ccx(ctrl, regN[i], regN_ctrl[i])
        c = Adder(c, regN_ctrl, reg2, ancil, reverse = False)
        #to restore regN_ctrl to 0
        for i in range(n):
            c.ccx(ctrl, regN[i], regN_ctrl[i])
        c = Adder(c, reg1, reg2, ancil, reverse=True)
        c.cx(ancil[0], ctrl)
        c = Adder(c, reg1, reg2, ancil, reverse=False)
        c.x(ctrl) #otherwise left to one
    else:
        c.x(ctrl) #otherwise left to one
        c = Adder(c, reg1, reg2, ancil, reverse=True)
        c.cx(ancil[0], ctrl)
        c = Adder(c, reg1, reg2, ancil, reverse=False)
        for i in range(n):
            i = n - 1 - i
            c.ccx(ctrl, regN[i], regN_ctrl[i])    
        c = Adder(c, regN_ctrl, reg2, ancil, reverse = True)
        for i in range(n):
            i = n - 1 - i
            c.ccx(ctrl, regN[i], regN_ctrl[i])
        c.cx(ancil[0], ctrl)
        c = Adder(c, regN, reg2, ancil, reverse = False)
        c = Adder(c, reg1, reg2, ancil, reverse = True)
    c.to_gate(label="Adder mod")
    return c

    
def Mult_mod(c, regx, reg1, ctrl, reg2, regN, regN_ctrl, ancil_ctrl, reverse=False):
    """
    Function: 
    It does ax mod n by doing \sum a x_{k}2^{k}.
    If ctrl is |1> and |x_{k}> is |1>, the register inputed in the Adder_mod is x_{k} as |x_{k} x_{k} ... x_{k} * n> 
    Args:
    c: circuit
    reg1: one of the number added of length n (reg1 > 0)
    reg2: the other number of length n, note that reg2 < N
    regN: register used to stored N used for the modulo
    regN_ctrl: register used to copy N
    ancil_mod: the n carrying qubits required and one extra controlled qubit 
    reverse (opt): bool to run the operation in reverse
    
    Return:
    The circuit which carries intrisically the qubit's state. 
    reg1 stays unchanged.
    reg2 becomes (reg1+reg2)%N for reverse = False 
    reg2 becomes reg2-reg1 if reg2-reg1>N else N+reg2-reg1 for reverse = True
    ancil_mod in fully set to 0
    reg_N remains unchanged
    reg_N_ctrl is set back to 0
    """    
    n = len(regx)
    if n != len(reg2):
        raise Exception("Currently only supports addition of equal bit length")
    try:
        ancil_mod = ancil_ctrl[0:n+2]
    except:
        raise Exception(f"Not enough ancilla with {len(ancil_ctrl)} given for {n+3} required")
    if not reverse:
        for i in range(0, n):
            ls = list(range(0, n))
            ls.pop(i - 1)
            for j in ls:
                c.ccx(ctrl, regx[i], reg1[j])
                c = Adder_mod(c, reg1, reg2, regN, regN_ctrl, ancil_mod)
                c.ccx(ctrl, regx[i], reg1[j])
        c.x(ctrl)
        for i in range(n):
            c.ccx(ctrl, regx[i], reg2[i])
        c.x(ctrl)        
    else:
        c.x(ctrl)
        #x is dealt with from 0 to n in indexing, but everything is reversed. but not here
        for i in range(n):
            i = n - 1 - i
            c.ccx(ctrl, regx[i], reg2[i])
        c.x(ctrl)
        for i in range(0, n):
            i = n - 1 - i
            ls = list(range(0, n))[::-1]
            ls.pop(i - 1)
            for j in ls:
                c.ccx(ctrl, regx[i], reg1[j])
                c = Adder_mod(c, reg1, reg2, regN, regN_ctrl, ancil_mod, reverse=True)
                c.ccx(ctrl, regx[i], reg1[j])  
    c.to_gate(label="Mult mod")
    return c
    

def Exp_mod(c, regx, reg1, regX, reg2, regN, regN_ctrl, ancil_ctrl):
    n = len(regX)
    c.x(regx)
    for i in range(n):
        c = Ctrl_mod(c, regx, reg1, regX[i], reg2, regN, regN_ctrl, ancil_ctrl, reverse=False)
        c = Ctrl_mod(c, reg1, regx, regX[i], reg2, regN, regN_ctrl, ancil_ctrl, reverse=True)
    c.to_gate(label="Exp mod")
    return c
    

    
    