import numpy as np
import qiskit as q
def nb_to_reg(nb_10, c, base=2, name="reg", size = None, classical = False):
    nb = np.base_repr(nb_10, base=base)
    if size is not None:
        nb = np.base_repr(nb_10, base=base, padding=size-len(nb))
    if classical:
        reg = q.ClassicalRegister(len(nb), name=name)
    else:
        reg = q.QuantumRegister(len(nb), name=name)
    c.add_register(reg)
    for i in range(len(nb)):
        if nb[i] == '1':
            c.x(reg[i])
    return reg, c

def keys_to_nb(keys, size, nb_keys = 2):
    ls = list(keys)[0].split(" ")
    if isinstance(size, int):
        size = nb_keys*[size]
    tmp = ls[0][::-1]
    ls = ls[1:]
    ls = [tmp[sum(size):]] + ls
    for i in range(nb_keys):
        s = size[i]
        i = nb_keys - i 
        ls = [tmp[s*(i-1):s*i]] + ls
        
    for i in range(nb_keys):
        ls[i] = int(ls[i], 2)
    return ls

def size(nb_list):
    ls = []
    if nb_list[-1] == 1:
        raise Exception('1 is not supported as trivial result')
    if not (nb_list[1] < nb_list[-1] or 0 <= nb_list[0]):
        raise Exception('The scheme restrains 0 <= a and b < N')
    n = len(np.base_repr(max(nb_list + [nb_list[0]+nb_list[1]])))
    return nb_list + [n]