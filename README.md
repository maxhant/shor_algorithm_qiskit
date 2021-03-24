# Shor's algorithm based on Vedral's version 
This algorithm is inteded to replicated the architecture in Vedral's paper (doi: 10.1103/PhysRevA.54.147).<br>
With Qiskit, it is coded up to the gate level using only C-NOT and Tofoli gates. It aims to using Qiskit's functions as little as possible. <br>
Therefore, it can serve as a basis for other quantum computer's API by simply replacing the two gates. <br>

<b>The algorithm is still under development as it cannot yet be run on the qasm simulator as it requires too much ressources.</b>

## Requirements
numpy <br>
qiskit <br>
pytest<br>
matplotlib<br>

## Installation
Clone <br>
Run pytest<br>
Use a jupyter script to try the functions<br>

## Visualization on Binder
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/maxhant/shor_algorithm_qiskit/main)
