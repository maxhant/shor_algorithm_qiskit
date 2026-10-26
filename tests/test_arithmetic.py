"""Exhaustive tests of the reversible arithmetic on computational basis states."""

from math import gcd

import pytest
from qiskit import QuantumCircuit, QuantumRegister

from shor_qiskit.arithmetic import (
    adder,
    adder_mod,
    bit_sum,
    carry,
    constant_bits,
    ctrl_mult_mod,
    mod_exp,
)

from .basis_sim import int_bits, load_int, read_int, run

N_BITS = 4
ALLOWED_GATES = {"x", "cx", "ccx"}


def only_allowed_gates(qc: QuantumCircuit) -> bool:
    return set(qc.count_ops()) <= ALLOWED_GATES


# --- building blocks ---------------------------------------------------------


@pytest.mark.parametrize("c_in", [0, 1])
@pytest.mark.parametrize("a", [0, 1])
@pytest.mark.parametrize("b", [0, 1])
def test_carry(c_in, a, b):
    r = QuantumRegister(4)
    qc = QuantumCircuit(r)
    load_int(qc, r[:3], c_in | a << 1 | b << 2)
    carry(qc, r[0], r[1], r[2], r[3])
    s = run(qc)
    assert (s[r[0]], s[r[1]]) == (c_in, a)
    assert s[r[2]] == a ^ b  # the sum bit is completed later by bit_sum
    assert s[r[3]] == int(a + b + c_in >= 2)


@pytest.mark.parametrize("c", [0, 1])
@pytest.mark.parametrize("a", [0, 1])
@pytest.mark.parametrize("b", [0, 1])
def test_bit_sum(c, a, b):
    r = QuantumRegister(3)
    qc = QuantumCircuit(r)
    load_int(qc, r, c | a << 1 | b << 2)
    bit_sum(qc, r[0], r[1], r[2])
    s = run(qc)
    assert (s[r[0]], s[r[1]], s[r[2]]) == (c, a, (a + b + c) % 2)


def adder_circuit(n):
    a = QuantumRegister(n, "a")
    b = QuantumRegister(n + 1, "b")
    c = QuantumRegister(n - 1, "c")
    return QuantumCircuit(a, b, c), a, b, c


@pytest.mark.parametrize("inverse", [False, True])
@pytest.mark.parametrize("n", [1, 2, 3, N_BITS])
def test_adder(n, inverse):
    for a_val in range(2**n):
        for b_val in range(2 ** (n + 1)):
            qc, a, b, c = adder_circuit(n)
            load_int(qc, a, a_val)
            load_int(qc, b, b_val)
            adder(qc, a, b, c, inverse=inverse)
            s = run(qc)
            expected = (b_val - a_val if inverse else a_val + b_val) % 2 ** (n + 1)
            assert read_int(s, b) == expected
            assert read_int(s, a) == a_val
            assert read_int(s, c) == 0


@pytest.mark.parametrize("inverse", [False, True])
def test_adder_with_constant_addend(inverse):
    n = N_BITS
    for a_val in range(2**n):
        for b_val in range(2**n):
            qc, _, b, c = adder_circuit(n)
            load_int(qc, b, b_val)
            adder(qc, constant_bits(a_val, n), b, c, inverse=inverse)
            expected = (b_val - a_val if inverse else a_val + b_val) % 2 ** (n + 1)
            assert read_int(run(qc), b) == expected


def mod_circuit(n):
    a = QuantumRegister(n, "a")
    b = QuantumRegister(n + 1, "b")
    c = QuantumRegister(n - 1, "c")
    t = QuantumRegister(1, "t")
    return QuantumCircuit(a, b, c, t), a, b, c, t[0]


@pytest.mark.parametrize("inverse", [False, True])
@pytest.mark.parametrize("N", range(2, 2**N_BITS))
def test_adder_mod(N, inverse):
    n = N_BITS
    for a_val in range(N):
        for b_val in range(N):
            qc, a, b, c, t = mod_circuit(n)
            load_int(qc, a, a_val)
            load_int(qc, b, b_val)
            adder_mod(qc, a, b, N, c, t, inverse=inverse)
            assert only_allowed_gates(qc)
            s = run(qc)
            expected = (b_val - a_val if inverse else a_val + b_val) % N
            assert read_int(s, b) == expected
            assert read_int(s, a) == a_val
            assert read_int(s, c) == 0 and s[t] == 0


def test_adder_mod_rejects_modulus_too_large():
    qc, a, b, c, t = mod_circuit(3)
    with pytest.raises(ValueError):
        adder_mod(qc, a, b, 8, c, t)


# --- multiplication and exponentiation ----------------------------------------


def mult_circuit(n):
    ctrl = QuantumRegister(1, "ctrl")
    x = QuantumRegister(n, "x")
    z = QuantumRegister(n + 1, "z")
    c = QuantumRegister(n - 1, "c")
    anc = QuantumRegister(2, "anc")
    return QuantumCircuit(ctrl, x, z, c, anc), ctrl[0], x, z, c, anc


@pytest.mark.parametrize("N", [3, 5, 7, 9, 11, 13, 15])
def test_ctrl_mult_mod(N):
    n = N.bit_length()
    for a_val in (a for a in range(1, N) if gcd(a, N) == 1):
        qc, ctrl, x, z, c, anc = mult_circuit(n)
        ctrl_mult_mod(qc, ctrl, x, z, a_val, N, c, anc[0], anc[1])
        assert only_allowed_gates(qc)
        for x_val in range(N):
            for ctrl_val in (0, 1):
                s = run(qc, {ctrl: ctrl_val, **int_bits(x, x_val)})
                expected = a_val * x_val % N if ctrl_val else x_val
                assert read_int(s, z) == expected
                assert read_int(s, x) == x_val
                assert read_int(s, c) == 0 and read_int(s, anc) == 0


@pytest.mark.parametrize("N", [7, 15])
def test_ctrl_mult_mod_inverse_undoes_forward(N):
    n = N.bit_length()
    a_val = 2
    for x_val in range(N):
        for ctrl_val in (0, 1):
            qc, ctrl, x, z, c, anc = mult_circuit(n)
            load_int(qc, [ctrl], ctrl_val)
            load_int(qc, x, x_val)
            ctrl_mult_mod(qc, ctrl, x, z, a_val, N, c, anc[0], anc[1])
            ctrl_mult_mod(qc, ctrl, x, z, a_val, N, c, anc[0], anc[1], inverse=True)
            s = run(qc)
            assert read_int(s, x) == x_val and read_int(s, z) == 0


@pytest.mark.parametrize("N", [15, 21])
def test_mod_exp(N):
    n = N.bit_length()
    n_exp = 4
    for a_val in (a for a in range(2, N) if gcd(a, N) == 1):
        e = QuantumRegister(n_exp, "e")
        x = QuantumRegister(n + 1, "x")
        z = QuantumRegister(n + 1, "z")
        c = QuantumRegister(n - 1, "c")
        anc = QuantumRegister(2, "anc")
        qc = QuantumCircuit(e, x, z, c, anc)
        result = mod_exp(qc, e, x, z, a_val, N, c, anc[0], anc[1])
        other = list(z) if result == list(x) else list(x)
        assert only_allowed_gates(qc)
        for e_val in range(2**n_exp):
            s = run(qc, int_bits(e, e_val))
            assert read_int(s, result) == pow(a_val, e_val, N)
            assert read_int(s, other) == 0
            assert read_int(s, e) == e_val
            assert read_int(s, c) == 0 and read_int(s, anc) == 0


def test_mod_exp_requires_coprime_base():
    e, x, z = QuantumRegister(2), QuantumRegister(5), QuantumRegister(5)
    c, anc = QuantumRegister(3), QuantumRegister(2)
    qc = QuantumCircuit(e, x, z, c, anc)
    with pytest.raises(ValueError):
        mod_exp(qc, e, x, z, 6, 15, c, anc[0], anc[1])
