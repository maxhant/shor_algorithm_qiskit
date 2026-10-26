"""Reversible modular arithmetic for Shor's algorithm.

This module implements the networks of

    V. Vedral, A. Barenco and A. Ekert,
    "Quantum networks for elementary arithmetic operations",
    Phys. Rev. A 54, 147 (1996), doi:10.1103/PhysRevA.54.147

using only X, CNOT and Toffoli gates, so it can be ported to another
quantum SDK by re-implementing three gates.

Conventions
-----------
* Registers are **little-endian**: index 0 is the least significant bit,
  the same convention Qiskit uses for measurement results.
* The classical numbers (the base ``a`` and the modulus ``N``) are not stored
  in quantum registers. Wherever the paper adds one of them from a register,
  each bit is *hard-wired* instead: a 0 bit drops its gates and a 1 bit uses a
  control qubit that is 1 exactly when the number should be added (or no
  control at all, for an unconditional addition). The adder only ever uses
  the addend as a control, so this produces the same unitary as the paper's
  networks while saving the "a", "N" and "temporary N" registers.
* Every function appends its gates to ``qc`` in place. Passing
  ``inverse=True`` appends the inverse network instead (for the adder that is
  subtraction). All gates used are self-inverse, so the inverse is the same
  gate list in reverse order.

Qubit budget for an ``n``-bit modulus (``N < 2**n``):

==============  =========  ==================================================
register        size       role
==============  =========  ==================================================
``x``           n + 1      multiplicand (the top qubit always stays 0)
``z``           n + 1      accumulator; top qubit is the overflow / sign bit
``carries``     n - 1      carry qubits of the ripple adder (all return to 0)
``t``           1          "a + b < N" flag of the modular adder
``flag``        1          ``ctrl AND x_k`` while adding ``a * 2**k mod N``
==============  =========  ==================================================

That is 3n + 3 work qubits, plus the counting register of Shor's algorithm.
"""

from __future__ import annotations

import functools
from collections.abc import Sequence
from typing import Union

from qiskit import QuantumCircuit
from qiskit.circuit import Qubit

#: How one bit of an addend is supplied: a qubit, or a classical constant.
Bit = Union[Qubit, bool]

__all__ = [
    "Bit",
    "carry",
    "bit_sum",
    "adder",
    "adder_mod",
    "ctrl_mult_mod",
    "mod_exp",
    "constant_bits",
    "work_qubits",
]


def _reversible(build):
    """Give a gate-building function an ``inverse`` keyword argument."""

    @functools.wraps(build)
    def wrapper(qc: QuantumCircuit, *args, inverse: bool = False, **kwargs):
        if not inverse:
            return build(qc, *args, **kwargs)
        block = qc.copy_empty_like()
        result = build(block, *args, **kwargs)
        qc.compose(block.inverse(), inplace=True)
        return result

    return wrapper


def _cx(qc: QuantumCircuit, control: Bit, target: Qubit) -> None:
    """CNOT whose control may be a classical constant."""
    if control is True:
        qc.x(target)
    elif control is not False:
        qc.cx(control, target)


def _ccx(qc: QuantumCircuit, control1: Bit, control2: Qubit, target: Qubit) -> None:
    """Toffoli whose first control may be a classical constant."""
    if control1 is True:
        qc.cx(control2, target)
    elif control1 is not False:
        qc.ccx(control1, control2, target)


def constant_bits(value: int, n: int, control: Bit = True) -> list[Bit]:
    """Hard-wire the classical number ``value`` as an ``n``-bit addend.

    Bits that are 1 become ``control`` (the constant ``True`` for an
    unconditional addition, or a qubit to make the addition controlled);
    bits that are 0 become ``False``.
    """
    if not 0 <= value < 2**n:
        raise ValueError(f"{value} does not fit in {n} bits")
    return [control if (value >> i) & 1 else False for i in range(n)]


@_reversible
def carry(qc: QuantumCircuit, c_in: Bit, a: Bit, b: Qubit, c_out: Qubit) -> None:
    """CARRY block: ``c_out ^= majority(c_in, a, b)`` and ``b ^= a``.

    The carry-in ``c_in`` and the addend bit ``a`` are only used as controls,
    so either may be a classical constant.
    """
    _ccx(qc, a, b, c_out)
    _cx(qc, a, b)
    _ccx(qc, c_in, b, c_out)


def bit_sum(qc: QuantumCircuit, c: Bit, a: Bit, b: Qubit) -> None:
    """SUM block: ``b ^= a ^ c``, i.e. ``b`` becomes the sum bit.

    The block is its own inverse.
    """
    _cx(qc, a, b)
    _cx(qc, c, b)


@_reversible
def adder(
    qc: QuantumCircuit,
    a: Sequence[Bit],
    b: Sequence[Qubit],
    carries: Sequence[Qubit],
) -> None:
    """Plain adder: ``b <- a + b`` on ``n + 1`` bits.

    Args:
        a: the ``n``-bit addend (qubits and/or constants), left unchanged.
        b: ``n + 1`` qubits; ``b[n]`` receives the carry out, so the result
            is ``(a + b) mod 2**(n+1)``.
        carries: ``n - 1`` qubits in state 0, returned to 0.

    With ``inverse=True`` this computes ``b <- (b - a) mod 2**(n+1)``; if
    ``b < 2**n`` on input, ``b[n]`` ends up as 1 exactly when ``a > b``.
    """
    n = len(a)
    if n < 1:
        raise ValueError("the addend must have at least one bit")
    if len(b) != n + 1:
        raise ValueError(f"b needs n + 1 = {n + 1} qubits, got {len(b)}")
    if len(carries) != n - 1:
        raise ValueError(f"the adder needs n - 1 = {n - 1} carry qubits, got {len(carries)}")

    # c[0] is the carry into the lowest bit, which is always 0.
    c: list[Bit] = [False, *carries, b[n]]

    for i in range(n):
        carry(qc, c[i], a[i], b[i], c[i + 1])
    _cx(qc, a[n - 1], b[n - 1])
    bit_sum(qc, c[n - 1], a[n - 1], b[n - 1])
    for i in reversed(range(n - 1)):
        carry(qc, c[i], a[i], b[i], c[i + 1], inverse=True)
        bit_sum(qc, c[i], a[i], b[i])


@_reversible
def adder_mod(
    qc: QuantumCircuit,
    a: Sequence[Bit],
    b: Sequence[Qubit],
    N: int,
    carries: Sequence[Qubit],
    t: Qubit,
) -> None:
    """Modular adder: ``b <- (a + b) mod N``.

    Requires ``0 <= a, b < N < 2**n``; the result is undefined otherwise.
    ``b`` has ``n + 1`` qubits and its top qubit must be 0. The ancilla ``t``
    and the carries start and end in 0. With ``inverse=True`` this computes
    ``b <- (b - a) mod N``.
    """
    n = len(a)
    if not 1 < N < 2**n:
        raise ValueError(f"N = {N} must satisfy 1 < N < 2**{n}")
    overflow = b[n]

    adder(qc, a, b, carries)                                   # b = a + b
    adder(qc, constant_bits(N, n), b, carries, inverse=True)   # b = a + b - N
    qc.cx(overflow, t)                                         # t = (a + b < N)
    adder(qc, constant_bits(N, n, control=t), b, carries)      # if t: b += N
    # b is now (a + b) mod N; t is still entangled with it. Uncompute t:
    adder(qc, a, b, carries, inverse=True)                     # b = result - a
    qc.x(overflow)                                             # negative iff t = 0
    qc.cx(overflow, t)
    qc.x(overflow)
    adder(qc, a, b, carries)                                   # b = (a + b) mod N


@_reversible
def ctrl_mult_mod(
    qc: QuantumCircuit,
    ctrl: Qubit,
    x: Sequence[Qubit],
    z: Sequence[Qubit],
    a: int,
    N: int,
    carries: Sequence[Qubit],
    t: Qubit,
    flag: Qubit,
) -> None:
    """Controlled modular multiplier.

    On input ``|ctrl, x, z=0>`` this produces::

        ctrl = 1:  |1, x, a * x mod N>
        ctrl = 0:  |0, x, x>

    using ``a * x = sum_k x_k * (a * 2**k mod N)``: for every bit ``x_k``
    the classical constant ``a * 2**k mod N`` is added with :func:`adder_mod`,
    conditioned on ``ctrl AND x_k``. Copying ``x`` when ``ctrl = 0`` is what
    lets :func:`mod_exp` swap registers unconditionally.

    Args:
        x: ``n`` qubits holding ``0 <= x < N``.
        z: ``n + 1`` qubits, initially 0.
        a: classical multiplier, with ``gcd(a, N) = 1`` for the operation
            to be reversible on the subspace that matters.
        t, flag: ancillas in state 0, returned to 0.
    """
    n = len(x)
    if len(z) != n + 1:
        raise ValueError(f"z needs n + 1 = {n + 1} qubits, got {len(z)}")

    for k in range(n):
        addend = (a << k) % N
        if addend == 0:
            continue
        qc.ccx(ctrl, x[k], flag)
        adder_mod(qc, constant_bits(addend, n, control=flag), z, N, carries, t)
        qc.ccx(ctrl, x[k], flag)

    # If ctrl = 0 the loop added nothing: copy x into z instead.
    qc.x(ctrl)
    for k in range(n):
        qc.ccx(ctrl, x[k], z[k])
    qc.x(ctrl)


def mod_exp(
    qc: QuantumCircuit,
    exponent: Sequence[Qubit],
    x: Sequence[Qubit],
    z: Sequence[Qubit],
    a: int,
    N: int,
    carries: Sequence[Qubit],
    t: Qubit,
    flag: Qubit,
) -> list[Qubit]:
    """Modular exponentiation: ``|e, 0> -> |e, a**e mod N>``.

    ``x`` and ``z`` are two ``(n + 1)``-qubit registers in state 0.
    ``x`` is set to 1, then for every exponent bit ``e_i``:

    1. ``z <- a_i * x mod N`` controlled by ``e_i``, with
       ``a_i = a**(2**i) mod N`` (copies ``x`` if ``e_i = 0``);
    2. the roles of ``x`` and ``z`` are swapped (a relabelling, no gates);
    3. the old register is cleared by running the multiplier backwards with
       the modular inverse ``a_i**-1 mod N``: ``x - a_i**-1 * (a_i * x) = 0``.

    Returns:
        The register (``x`` or ``z``) that ends up holding ``a**e mod N``;
        the other one is back to 0.
    """
    from math import gcd

    if gcd(a, N) != 1:
        raise ValueError(f"a = {a} and N = {N} must be coprime")
    n = len(x) - 1
    x, z = list(x), list(z)
    qc.x(x[0])
    for i, e in enumerate(exponent):
        a_i = pow(a, 2**i, N)
        if a_i == 1:
            # Multiplying by 1 is the identity: skip the whole block.
            continue
        ctrl_mult_mod(qc, e, x[:n], z, a_i, N, carries, t, flag)
        x, z = z, x
        ctrl_mult_mod(qc, e, x[:n], z, pow(a_i, -1, N), N, carries, t, flag, inverse=True)
    return x


def work_qubits(N: int) -> int:
    """Number of work qubits (excluding the counting register) used for ``N``."""
    n = N.bit_length()
    return 3 * n + 3
