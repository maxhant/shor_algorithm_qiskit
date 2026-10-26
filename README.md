# Shor's algorithm based on Vedral's arithmetic

An implementation of Shor's factoring algorithm in Qiskit, with the modular
exponentiation built gate by gate from the reversible networks of

> V. Vedral, A. Barenco and A. Ekert, *Quantum networks for elementary arithmetic
> operations*, Phys. Rev. A **54**, 147 (1996), [doi:10.1103/PhysRevA.54.147](https://doi.org/10.1103/PhysRevA.54.147)

All the arithmetic uses only X, CNOT and Toffoli gates and none of Qiskit's
library circuits, so it can be ported to another quantum SDK by
re-implementing those three gates. The quantum Fourier transform is likewise
written out with Hadamard, controlled-phase and SWAP gates.

It factors 15 and 21 on a laptop in a few seconds.

## Quick start

```bash
git clone https://github.com/maxhant/shor_algorithm_qiskit.git
cd shor_algorithm_qiskit
pip install -e ".[test,notebook]"
pytest
```

```python
from shor_qiskit import factor, find_order

factor(21, seed=3, allow_lucky_gcd=False)
# ShorResult(N=21, factors=(3, 7), method='order finding', a=19, order=6, attempts=2)

find_order(7, 15)[0]
# 4
```

The notebook [`Project_visualization.ipynb`](Project_visualization.ipynb) walks
through every block, from a one-bit SUM to factoring, with circuits you can
modify and run.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/maxhant/shor_algorithm_qiskit/HEAD?labpath=Project_visualization.ipynb)

## How it is built

| Block | Function | What it does |
|---|---|---|
| SUM | `bit_sum` | `b ^= a ^ c`, the sum bit of a one-bit addition |
| CARRY | `carry` | `c_out ^= majority(c_in, a, b)` |
| ADDER | `adder` | `b <- a + b` (or `b - a` with `inverse=True`) |
| ADDER MOD N | `adder_mod` | `b <- (a + b) mod N` for `0 <= a, b < N` |
| C-MULT MOD N | `ctrl_mult_mod` | `z <- a * x mod N` if the control is 1, `z <- x` otherwise |
| MODULAR EXP | `mod_exp` | `|e>|0> -> |e>|a^e mod N>` |
| QFT | `qft`, `inverse_qft` | quantum Fourier transform on the counting register |
| Shor | `order_finding_circuit`, `find_order`, `factor` | order finding, continued fractions and the classical reduction to factoring |

Every arithmetic function appends to a circuit in place and takes
`inverse=True` to append the inverse network. Registers are little-endian
(qubit 0 is the least significant bit), as in the rest of Qiskit.

### Differences from the paper

The block structure is the paper's, with one resource optimisation: the
classical numbers `a` and `N` are not stored in quantum registers. The adder
only ever uses its addend as a control, so each bit of a classical addend is
hard-wired instead: a 0 bit drops its gates, and a 1 bit uses either no
control (to add `N`) or a single control qubit (to add `N` only when the
modular adder's flag is set, or to add `a * 2^k mod N` only when
`ctrl AND x_k`). This removes the paper's `a`, `N` and temporary-`N`
registers and the carry qubit that is always 0.

For an `n`-bit `N` the circuit uses `3n + 3` work qubits plus a `2n`-qubit
counting register: 23 qubits for N = 15 and 28 for N = 21. Those are simulated
with Aer's matrix-product-state method, which is fast here because the state
stays weakly entangled (its Schmidt rank is at most the order `r`).

## Tests

`pytest` checks each arithmetic block exhaustively on all basis-state inputs
(with a small classical simulator, since X/CNOT/Toffoli circuits only permute
basis states), the QFT against Qiskit's `QFTGate`, the exact output
distribution of order finding for N = 15, and end-to-end factoring of 15
and 21. The suite runs in about 30 seconds.
