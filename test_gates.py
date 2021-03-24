import qiskit as q
import numpy as np
from gates import *
from handler import nb_to_reg, keys_to_nb, size
import pytest
from qiskit import Aer
backend = Aer.get_backend('qasm_simulator')

def test_Carry():
    c = q.QuantumCircuit()
    qubits = q.QuantumRegister(4, name='qubits')
    c.add_register(qubits)
    c1 = np.random.randint(2)
    n1 = np.random.randint(2)
    n2 = np.random.randint(2)
    if c1 == 1:
        c.x(qubits[0])
    if n1 == 1:
        c.x(qubits[1])
    if n2 == 1:
        c.x(qubits[2])
    if (n1+n2+c1)>1:
        carr = 1
    else:
        carr = 0
    c = Carry(c, qubits)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)
    out = list(out)[0]
    assert (n1+n2+c1)%2 == int(out[1])
    assert carr == int(out[0])
    assert c1 == int(out[3])
    assert n1 == int(out[2])

def test_Sum():
    c = q.QuantumCircuit()
    qubits = q.QuantumRegister(3, name='qubits')
    c.add_register(qubits)
    c1 = np.random.randint(2)
    n1 = np.random.randint(2)
    if c1 == 1:
        c.x(qubits[0])
    if n1 == 1:
        c.x(qubits[1])
    c = Sum(c, qubits)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)
    out = list(out)[0]
    assert c1 + n1 == int(out[0])
    assert c1 == int(out[2])
    assert n1 == int(out[1])
    

def test_Adder():
    n = 4
    nb_max = int('1'*n, 2)
    n1 = np.random.randint(1, int(nb_max/2))
    n2 = np.random.randint(1, int(nb_max/2))
    c = q.QuantumCircuit()
    reg1, c = nb_to_reg(n1, c, name="a_{reg}", size=n)
    reg2, c = nb_to_reg(n2, c, name="b_{reg}", size=n)
    ancil = q.QuantumRegister(n+1, name='ancillae')
    meas1 = q.ClassicalRegister(n, name='a_{meas}')
    meas2 = q.ClassicalRegister(n, name='b_{meas}')
    c.add_register(ancil, meas1, meas2)

    c = Adder(c, reg1, reg2, ancil)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)
    vals = keys_to_nb(out, size=n, nb_keys = 2)
    assert n1 + n2 == vals[1]
    assert n1 == vals[0]
    assert n*'0' == vals[2][1:]
    if n1+n2>=nb_max:
        assert '1' == vals[2][0]
    else:
        assert '0' == vals[2][0]


def test_Adder_r():
    n = 4
    nb_max = int('1'*n, 2)
    n1, n2 = np.sort([np.random.randint(1, int(nb_max/2)), np.random.randint(1, int(nb_max/2))])
    c = q.QuantumCircuit()
    reg1, c = nb_to_reg(n1, c, name="a_{reg}", size=n)
    reg2, c = nb_to_reg(n2, c, name="b_{reg}", size=n)
    ancil = q.QuantumRegister(n+1, name='ancillae')
    meas1 = q.ClassicalRegister(n, name='a_{meas}')
    meas2 = q.ClassicalRegister(n, name='b_{meas}')
    c.add_register(ancil, meas1, meas2)

    c = Adder(c, reg1, reg2, ancil, reverse=True)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)           
    vals = keys_to_nb(out, size=n, nb_keys=2)
    assert n2 - n1 == vals[1]
    assert n1 == vals[0]
    assert n*'0' == vals[2][1:] 
    
def test_Adder_mod():
    n = 4
    nb_max = int('1'*n, 2)
    n1, N, n = size([np.random.randint(1, int(nb_max/2)), np.random.randint(2, int(nb_max))])
    n2 = np.random.randint(1, N)
    c = q.QuantumCircuit()  
    reg1, c = nb_to_reg(n1, c, name="a_{reg}", size=n)
    reg2, c = nb_to_reg(n2, c, name="b_{reg}", size=n)
    regN, c = nb_to_reg(N, c, name="capn_{reg}", size=n)
    regN_ctrl = q.QuantumRegister(n, name='n_ctrl')
    ancil = q.QuantumRegister(n+2, name='ancillae')
    meas1 = q.ClassicalRegister(n, name='a_{meas}')
    meas2 = q.ClassicalRegister(n, name='b_{meas}')
    c.add_register(regN_ctrl, ancil, meas1, meas2)

    c = Adder_mod(c, reg1, reg2, regN, regN_ctrl, ancil)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)
    vals = keys_to_nb(out, size=n, nb_keys = 4)
    assert (n2+n1)%N == vals[1]
    assert n1 == vals[0]
    assert N == vals[2]
    assert 0 == vals[3]
    assert (n+2)*'0' == vals[4]
    
def test_Adder_mod_r():
    n = 4
    nb_max = int('1'*n, 2)
    n1, N, n = size([np.random.randint(1, int(nb_max/2)), np.random.randint(2, int(nb_max))])
    n2 = np.random.randint(1, N)
    c = q.QuantumCircuit()  
    reg1, c = nb_to_reg(n1, c, name="a_{reg}", size=n)
    reg2, c = nb_to_reg(n2, c, name="b_{reg}", size=n)
    regN, c = nb_to_reg(N, c, name="capn_{reg}", size=n)
    regN_ctrl = q.QuantumRegister(n, name='n_ctrl')
    ancil = q.QuantumRegister(n+2, name='ancillae')
    meas1 = q.ClassicalRegister(n, name='a_{meas}')
    meas2 = q.ClassicalRegister(n, name='b_{meas}')
    c.add_register(regN_ctrl, ancil, meas1, meas2)

    c = Adder_mod(c, reg1, reg2, regN, regN_ctrl, ancil, reverse=True)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)
    vals = keys_to_nb(out, size=n, nb_keys=4)
    if n2 - n1 != (n2-n1)%N:
        assert n2-n1+N == vals[1] 
    elif n2 - n1 < 0:
        assert n1 - n2 == vals[1]
    else:
        assert n2 - n1 == vals[1]
    assert n1 == vals[0]
    assert N == vals[2]
    assert 0 == vals[3]
    assert (n+2)*'0' == vals[4]
    
def test_Mult_mod():
    backend = Aer.get_backend('qasm_simulator')
    c = q.QuantumCircuit()
    """
    x is the exponent and it is initialized with H gate, thus it covers all values from 0 to 2^{n} - 1
    a is the guess done
    N is the coprime number
    """
    x, N, n = size([3, 7])
    regx, c = nb_to_reg(x, c, name="a_{reg}", size=n) #I dont think a is inputed here
    reg2 = q.QuantumRegister(n, name='reg_0')
    c.add_register(reg2)
    reg1 = q.QuantumRegister(n, name='reg_1')
    regN, c = nb_to_reg(N, c, name="capn_{reg}", size=n)
    regN_ctrl = q.QuantumRegister(n, name='n_ctrl')
    regX, c = nb_to_reg(1, c, name="x_{reg}", size=1)
    ancil = q.QuantumRegister(n+2, name='ancillae')
    meas1 = q.ClassicalRegister(n, name='a_{meas}')
    meas2 = q.ClassicalRegister(n, name='b_{meas}')
    c.add_register(reg1, regN_ctrl, ancil, meas1, meas2)

    c = Mult_mod(c, regx, reg1, regX[0], reg2, regN, regN_ctrl, ancil, reverse=False)
    c.measure_all()
    job = q.execute(c, backend=backend, shots=1)
    out = job.result().get_counts(c)
    vals = keys_to_nb(out, size=n, nb_keys=4)
    assert True
