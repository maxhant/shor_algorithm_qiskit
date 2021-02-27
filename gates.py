import qiskit as q
import numpy as np

def Carry(c, idx, reverse = False):
    if not len(idx) == 4:
        raise Exception("Carry gate takes 4 qubits")
    if reverse:
        c.ccx(idx[0],idx[2],idx[3])
        c.cx(idx[1],idx[2])
        c.ccx(idx[1],idx[2],idx[3])
    else:
        c.ccx(idx[1],idx[2],idx[3])
        c.cx(idx[1],idx[2])
        c.ccx(idx[0],idx[2],idx[3])
#     c.barrier()
    return c

def Sum(c, idx, reverse = False):
    if not len(idx) == 3:
        raise Exception("Carry gate takes 4 qubits")
    if reverse:
        c.cx(idx[0], idx[2])
        c.cx(idx[1], idx[2])        
    else:
        c.cx(idx[1], idx[2])
        c.cx(idx[0], idx[2])
#     c.barrier()
    return c

def Adder(c, reg1, reg2, ancil_adder, reverse = False):
    n = len(reg1)
    try:
        ancil = ancil_adder[0:n+1]
    except: 
        raise Exception(f" {len(ancil_adder)} ancillaes when {n+1} required.")
    if n != len(reg2):
        raise Exception("Currently only supports addition of equal bit length")        
    if reverse:
        for i in range(1, n):
            i = n - i
            c = Sum(c, [ancil[i+1], reg1[i], reg2[i]], reverse=True)            
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]], reverse=False)
        c = Sum(c, [ancil[1], reg1[0], reg2[0]], reverse=True)
        c.cx(reg1[0], reg2[0])
        for i in range(n):
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]], reverse=True)          
    else:
        for i in range(n):
            i = (n-1) - i
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]])
        c.cx(reg1[0], reg2[0])
        c = Sum(c, [ancil[1], reg1[0], reg2[0]])
        for i in range(1, n):
            c = Carry(c, [ancil[i+1], reg1[i], reg2[i], ancil[i]], reverse=True)
            c = Sum(c, [ancil[i+1], reg1[i], reg2[i]])
    return c

def Adder_mod(c, reg1, reg2, regN, regN_ctrl, ancil_mod):
    n = len(reg1)
    ancil = ancil_mod[0:n+1]
    if n != len(reg2):
        raise Exception("Currently only supports addition of equal bit length")
    try:
        ctrl = ancil_mod[n+1]
    except:
        raise Exception(f"Not enough ancilla with {len(ancil_mod)} given for {n+2} required")
    c = Adder(c, reg1, reg2, ancil, reverse = False)
    c.barrier()
    c = Adder(c, regN, reg2, ancil, reverse = True)
    c.barrier()
    c.cx(ancil[0], ctrl)
    c.barrier()
    for i in range(n):
        c.ccx(ctrl, regN[i], regN_ctrl[i])
    c.barrier()
    c = Adder(c, regN_ctrl, reg2, ancil, reverse = False)
    c.barrier() 
    #to restore regN_ctrl to 0
    for i in range(n):
        c.ccx(ctrl, regN[i], regN_ctrl[i])
    c = Adder(c, reg1, reg2, ancil, reverse=True)
    c.cx(ancil[0], ctrl)
    c = Adder(c, reg1, reg2, ancil, reverse=False)
    c.x(ctrl) #otherwise left to one

    return c

def Ctrl_zero(c, ctrl, regN, regN_ctrl):
    n = len(regN)
    if len(regN_ctrl) != n:
        raise Exception(f'Need {n} ancillae to copy N')
    # Ancils are all ones, then they are flipped if regN[i] is 1 so they become 
    # exact opposite
    for i in range(n):
        c.ccx(ctrl, regN[i], regN_ctrl[i])
    return c
    
    

    
    

    
    