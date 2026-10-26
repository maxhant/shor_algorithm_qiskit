"""Tests of the QFT, the order-finding circuit and the full factoring routine."""

from fractions import Fraction

import numpy as np
import pytest
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFTGate
from qiskit.quantum_info import Operator
from qiskit_aer import AerSimulator

from shor_qiskit import factor, find_order, order_finding_circuit, order_from_measurement
from shor_qiskit.qft import inverse_qft, qft


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_qft_matches_qiskit(n):
    qc = QuantumCircuit(n)
    qft(qc, qc.qubits)
    assert Operator(qc).equiv(Operator(QFTGate(n)))


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_inverse_qft_is_inverse(n):
    qc = QuantumCircuit(n)
    qft(qc, qc.qubits)
    inverse_qft(qc, qc.qubits)
    assert Operator(qc).equiv(Operator(np.eye(2**n)))


def test_order_finding_distribution_is_exact_for_15():
    """For a = 7, N = 15 (order 4), 3 counting qubits must give 0, 2, 4, 6 equally.

    Uses an exact statevector, independently of the default MPS simulator.
    """
    qc = order_finding_circuit(7, 15, n_count=3).remove_final_measurements(inplace=False)
    qc.save_probabilities(qubits=range(3))
    simulator = AerSimulator(method="statevector")
    probs = simulator.run(transpile(qc, simulator)).result().data()["probabilities"]
    peaks = {y: p for y, p in enumerate(probs) if p > 1e-9}
    assert set(peaks) == {0, 2, 4, 6}
    assert all(p == pytest.approx(0.25) for p in peaks.values())


def test_order_from_measurement():
    # r = 6, n_count = 10: y close to s * 1024 / 6
    assert order_from_measurement(round(1024 / 6), 10, 2, 21) == 6
    # s = 2 gives 1/3, a divisor of the order, and the multiples recover 6
    assert order_from_measurement(round(2 * 1024 / 6), 10, 2, 21) == 6
    assert order_from_measurement(0, 10, 2, 21) is None
    assert Fraction(683, 1024).limit_denominator(21) == Fraction(2, 3)


@pytest.mark.parametrize(("a", "N", "order"), [(2, 15, 4), (7, 15, 4), (4, 15, 2), (2, 21, 6), (5, 21, 6)])
def test_find_order(a, N, order):
    r, counts = find_order(a, N, shots=64, seed=1)
    assert r == order
    assert sum(counts.values()) == 64


def test_order_finding_rejects_bad_base():
    with pytest.raises(ValueError):
        order_finding_circuit(6, 15)


@pytest.mark.parametrize(("N", "factors"), [(15, (3, 5)), (21, (3, 7))])
def test_factor_runs_the_quantum_part(N, factors):
    result = factor(N, seed=3, allow_lucky_gcd=False)
    assert result.success
    assert result.factors == factors
    assert result.method == "order finding"
    assert pow(result.a, result.order, N) == 1


def test_factor_classical_shortcuts():
    assert factor(22).method == "even"
    assert factor(9).factors == (3, 3)
    assert factor(15, a=6).method == "gcd(a, N)"
    with pytest.raises(ValueError):
        factor(13)
