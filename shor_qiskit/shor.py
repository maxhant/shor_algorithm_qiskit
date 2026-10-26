"""Shor's algorithm: quantum order finding plus the classical reduction to factoring."""

from __future__ import annotations

import random
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile

from .arithmetic import mod_exp
from .qft import inverse_qft

__all__ = [
    "order_finding_circuit",
    "order_from_measurement",
    "find_order",
    "factor",
    "ShorResult",
    "default_simulator",
]


def order_finding_circuit(a: int, N: int, n_count: int | None = None) -> QuantumCircuit:
    """Build the order-finding circuit for ``a`` modulo ``N``.

    Layout: a counting register of ``n_count`` qubits (default ``2n`` for an
    ``n``-bit ``N``) put in uniform superposition, the modular exponentiation
    ``|e>|0> -> |e>|a**e mod N>``, an inverse QFT on the counting register and
    a measurement of it. The measured integer ``y`` is close to ``s * 2**n_count / r``
    for a random ``s``, where ``r`` is the order of ``a``.
    """
    if N < 3:
        raise ValueError("N must be at least 3")
    if not 1 < a < N or gcd(a, N) != 1:
        raise ValueError(f"a must satisfy 1 < a < N and gcd(a, N) = 1, got a = {a}")
    n = N.bit_length()
    n_count = 2 * n if n_count is None else n_count

    count = QuantumRegister(n_count, "count")
    x = QuantumRegister(n + 1, "x")
    z = QuantumRegister(n + 1, "z")
    carries = QuantumRegister(n - 1, "carry")
    anc = QuantumRegister(2, "anc")
    result = ClassicalRegister(n_count, "y")
    qc = QuantumCircuit(count, x, z, carries, anc, result, name=f"order_finding({a}, {N})")

    qc.h(count)
    mod_exp(qc, count, x, z, a, N, carries, anc[0], anc[1])
    inverse_qft(qc, count)
    qc.measure(count, result)
    return qc


def order_from_measurement(y: int, n_count: int, a: int, N: int) -> int | None:
    """Recover the order of ``a`` mod ``N`` from one measured value ``y``.

    The continued-fraction expansion of ``y / 2**n_count`` gives ``s / r`` in
    lowest terms, i.e. a divisor of the order when ``gcd(s, r) > 1``, so small
    multiples of the denominator are also tried. Returns ``None`` on failure.
    """
    if y == 0:
        return None
    denominator = Fraction(y, 2**n_count).limit_denominator(N).denominator
    multiple = denominator
    while multiple < N:
        if pow(a, multiple, N) == 1:
            return multiple
        multiple += denominator
    return None


def default_simulator():
    """Aer's matrix-product-state simulator.

    The state after the modular exponentiation has low entanglement (the
    Schmidt rank is at most the order ``r``), so this method handles the
    23 qubits needed for N = 15 and 28 for N = 21 in seconds, where a
    statevector would need gigabytes.
    """
    from qiskit_aer import AerSimulator

    return AerSimulator(method="matrix_product_state")


def find_order(
    a: int,
    N: int,
    *,
    shots: int = 32,
    n_count: int | None = None,
    simulator=None,
    seed: int | None = None,
) -> tuple[int | None, dict[str, int]]:
    """Run the order-finding circuit and return ``(order or None, counts)``."""
    simulator = simulator or default_simulator()
    qc = order_finding_circuit(a, N, n_count)
    n_count = qc.num_clbits
    compiled = transpile(qc, simulator, optimization_level=0)
    counts = simulator.run(compiled, shots=shots, seed_simulator=seed).result().get_counts()
    orders = []
    for bitstring in sorted(counts, key=counts.get, reverse=True):
        r = order_from_measurement(int(bitstring, 2), n_count, a, N)
        if r is not None:
            orders.append(r)
    # Every valid candidate is a multiple of the true order; the smallest one wins.
    return (min(orders) if orders else None), counts


@dataclass
class ShorResult:
    """Outcome of :func:`factor`."""

    N: int
    factors: tuple[int, int] | None
    method: str
    a: int | None = None
    order: int | None = None
    counts: Mapping[str, int] = field(default_factory=dict, repr=False)
    attempts: int = 0

    @property
    def success(self) -> bool:
        return self.factors is not None


def _is_probable_prime(N: int) -> bool:
    if N < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if N % p == 0:
            return N == p
    d, s = N - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(b, d, N)
        if x in (1, N - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, N)
            if x == N - 1:
                break
        else:
            return False
    return True


def _perfect_power(N: int) -> tuple[int, int] | None:
    for k in range(2, N.bit_length() + 1):
        root = round(N ** (1 / k))
        for b in (root - 1, root, root + 1):
            if b > 1 and b**k == N:
                return b, k
    return None


def factor(
    N: int,
    a: int | None = None,
    *,
    shots: int = 32,
    max_attempts: int = 10,
    n_count: int | None = None,
    simulator=None,
    seed: int | None = None,
    allow_lucky_gcd: bool = True,
) -> ShorResult:
    """Find a non-trivial factor of ``N`` with Shor's algorithm.

    Classical pre-checks handle the cases order finding cannot: even ``N``,
    primes (rejected) and perfect powers. Then, for a random ``a`` (or the
    one given), the order ``r`` is found on the simulator and
    ``gcd(a**(r/2) +- 1, N)`` is tried. This needs ``r`` even and
    ``a**(r/2) != -1 mod N``; otherwise another ``a`` is drawn.

    Args:
        a: base to use; if given, only this base is tried.
        allow_lucky_gcd: when a random ``a`` shares a factor with ``N``,
            return it directly (``True``) or draw another ``a`` so that the
            quantum part always runs (``False``).
        seed: makes the choice of ``a`` and the simulation reproducible.
    """
    if N < 4:
        raise ValueError("N must be a composite number >= 4")
    if N % 2 == 0:
        return ShorResult(N, (2, N // 2), "even")
    if _is_probable_prime(N):
        raise ValueError(f"{N} is prime")
    power = _perfect_power(N)
    if power is not None:
        b, _ = power
        return ShorResult(N, (b, N // b), "perfect power")

    rng = random.Random(seed)
    bases = [a] if a is not None else None
    last = ShorResult(N, None, "order finding failed")
    for attempt in range(1, max_attempts + 1):
        if bases is not None:
            if attempt > len(bases):
                break
            base = bases[attempt - 1]
        else:
            base = rng.randrange(2, N - 1)
        g = gcd(base, N)
        if g > 1:
            if allow_lucky_gcd or a is not None:
                return ShorResult(N, (g, N // g), "gcd(a, N)", a=base, attempts=attempt)
            continue

        r, counts = find_order(
            base, N, shots=shots, n_count=n_count, simulator=simulator,
            seed=None if seed is None else seed + attempt,
        )
        last = ShorResult(N, None, "order finding failed", base, r, counts, attempt)
        if r is None or r % 2:
            continue
        half = pow(base, r // 2, N)
        if half == N - 1:
            continue
        for candidate in (gcd(half - 1, N), gcd(half + 1, N)):
            if 1 < candidate < N:
                p, q = sorted((candidate, N // candidate))
                return ShorResult(N, (p, q), "order finding", base, r, counts, attempt)
    return last
