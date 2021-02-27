import qiskit as q
import numpy as np
from gates import *
from handler import nb_to_reg, keys_to_nb
import pytest
from qiskit import Aer
backend = Aer.get_backend('qasm_simulator')

def test_Adder():
    for N in range(1, 8):
        nb_max = int('1'*N, 2)
        for n1 in range(1,int(nb_max/2)):
            for n2 in range(1,int(nb_max/2)):
                c = q.QuantumCircuit()
                reg1, c = nb_to_reg(n1, c, name="a_{reg}", size=N)
                reg2, c = nb_to_reg(n2, c, name="b_{reg}", size=N)
                ancil = q.QuantumRegister(N+1, name='ancillae')
                meas1 = q.ClassicalRegister(N, name='a_{meas}')
                meas2 = q.ClassicalRegister(N, name='b_{meas}')
                c.add_register(ancil, meas1, meas2)

                c = Adder(c, reg1, reg2, ancil)
                c.measure_all()
                job = q.execute(c, backend=backend, shots=1)
                out = job.result().get_counts(c)           
                assert n1 + n2 == keys_to_nb(out, size=N)[1]