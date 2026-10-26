"""Quantum Fourier transform written out gate by gate (H, controlled-phase, SWAP)."""

from __future__ import annotations

from collections.abc import Sequence
from math import pi

from qiskit import QuantumCircuit
from qiskit.circuit import Qubit

__all__ = ["qft", "inverse_qft"]


def qft(qc: QuantumCircuit, qubits: Sequence[Qubit]) -> None:
    """Append the QFT on ``qubits`` (little-endian, same convention as Qiskit's ``QFTGate``)."""
    n = len(qubits)
    for j in reversed(range(n)):
        qc.h(qubits[j])
        for m in reversed(range(j)):
            qc.cp(pi / 2 ** (j - m), qubits[m], qubits[j])
    for i in range(n // 2):
        qc.swap(qubits[i], qubits[n - 1 - i])


def inverse_qft(qc: QuantumCircuit, qubits: Sequence[Qubit]) -> None:
    """Append the inverse QFT on ``qubits``."""
    n = len(qubits)
    for i in range(n // 2):
        qc.swap(qubits[i], qubits[n - 1 - i])
    for j in range(n):
        for m in range(j):
            qc.cp(-pi / 2 ** (j - m), qubits[m], qubits[j])
        qc.h(qubits[j])
