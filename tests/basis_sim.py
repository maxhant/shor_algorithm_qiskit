"""A tiny simulator for circuits made of X, CNOT, Toffoli and SWAP.

Such circuits map computational basis states to basis states, so they can be
checked exhaustively by tracking one bit per qubit instead of 2**n amplitudes.
"""

from __future__ import annotations

from collections.abc import Sequence

from qiskit import QuantumCircuit
from qiskit.circuit import Qubit


def load_int(qc: QuantumCircuit, qubits: Sequence[Qubit], value: int) -> None:
    """Prepare ``|value>`` (little-endian) on qubits that start in |0>."""
    if value >= 2 ** len(qubits):
        raise ValueError(f"{value} does not fit in {len(qubits)} qubits")
    for i, q in enumerate(qubits):
        if (value >> i) & 1:
            qc.x(q)


def run(qc: QuantumCircuit, initial: dict[Qubit, int] | None = None) -> dict[Qubit, int]:
    """Run ``qc`` and return the final bit of every qubit.

    Every qubit starts in 0 unless ``initial`` gives another bit for it.
    """
    state = {q: 0 for q in qc.qubits}
    state.update(initial or {})
    for instruction in qc.data:
        name, qs = instruction.operation.name, instruction.qubits
        if name == "barrier":
            continue
        if name == "x":
            state[qs[0]] ^= 1
        elif name == "cx":
            state[qs[1]] ^= state[qs[0]]
        elif name == "ccx":
            state[qs[2]] ^= state[qs[0]] & state[qs[1]]
        elif name == "swap":
            state[qs[0]], state[qs[1]] = state[qs[1]], state[qs[0]]
        else:
            raise ValueError(f"gate {name!r} is not classical-reversible")
    return state


def int_bits(qubits: Sequence[Qubit], value: int) -> dict[Qubit, int]:
    """Initial bits for :func:`run` that hold ``value`` (little-endian) on ``qubits``."""
    return {q: (value >> i) & 1 for i, q in enumerate(qubits)}


def read_int(state: dict[Qubit, int], qubits: Sequence[Qubit]) -> int:
    """Read a little-endian integer from a state returned by :func:`run`."""
    return sum(state[q] << i for i, q in enumerate(qubits))
