# Candidate exact characterization of singular-descriptor port-Hamiltonian realizations

Not peer reviewed. Not verified by a human expert.

Maintainer: Bruno L Yamamoto

This repository contains a candidate exact characterization and construction for realizing regular descriptor systems in port-Hamiltonian form.

This proof was found by GPT Astra Pro. It has not been checked by a human expert yet, and the result may already be known.


## Origin of the problem

This repository addresses the singular-descriptor extension of modern constructive port-Hamiltonian realization theory.

Relevant sources include:

- C. Beattie, V. Mehrmann, and H. Xu, *Port-Hamiltonian realizations of non-minimal linear time-invariant systems*, Mathematics of Control, Signals, and Systems 38 (2026), 1–41. DOI: `10.1007/s00498-025-00432-w`. This paper treats nonminimal standard LTI systems and identifies extension to descriptor systems with singular $E$ as an open direction.
- K. Cherifi, H. Gernandt, and D. Hinsen, *The difference between port-Hamiltonian, passive and positive real descriptor systems*, Mathematics of Control, Signals, and Systems 36 (2024), 451–482. DOI: `10.1007/s00498-023-00373-2`.
- D. Chu and V. Mehrmann, *Port-Hamiltonian Realizations of Positive Real Descriptor Systems*, arXiv:2408.14115 (2024), giving strong realization results under controllability/observability assumptions.

## Main result

For a regular descriptor realization

$$
E\dot x=Ax+Bu,\qquad y=Cx+Du,
$$

the proof claims a finite exact characterization of port-Hamiltonian realizability under the equivalence conventions stated in `PROOF.md`.

In particular, it claims:

- descriptor index greater than two is an obstruction;
- for index at most two, realizability reduces to affine KYP/LMI feasibility plus a finite determinant condition on the residual algebraic block;
- the determinant condition is decidable by a finite affine-hull/interpolation argument;
- nonminimal realizations, singular symmetric feedthrough, and singular energy multipliers are retained;
- the resulting matrices can be completed to the standard port-Hamiltonian dissipation form.

The full statement and proof are in `PROOF.md`.

## GitHub rendering

The Markdown uses GitHub-native MathJax syntax: `$...$` for inline mathematics and `$$...$$` for display mathematics.

## Repository contents

- `PROOF.md` — theorem statements, candidate proof, examples, literature discussion, and verification checklist.
- `code/` — reference implementation and verification scripts.
- `checks/` — exact certificates and saved solver inputs/outputs used by the checks.
- `requirements.txt` — Python dependency for the verification scripts.

## Reproducing the checks

Install the Python dependency:

```sh
python -m pip install -r requirements.txt
```

The independent certificate checker requires only SymPy:

```sh
python code/verify_certificates.py
```

A larger exact regression suite is also included:

```sh
python code/run_tests.py
```

`run_tests.py` additionally requires a system installation of the Z3 shared library (`libz3`). The script searches for it automatically; if needed, set `Z3_LIBRARY` to the full path of the shared library.

## Review and novelty status

This candidate proof has not yet undergone independent peer review or formal verification. Corrections, if needed, should be made transparently through the repository history.

The basic descriptor KYP-plus-kernel characterization is prior art. The candidate contribution here is the finite block reduction, determinant criterion, and constructive characterization stated in `PROOF.md`. A literature search did not find this exact theorem, but this is evidence of apparent novelty, not a certified priority claim.

## Scope

The claims concern only the system class, equivalence notions, and port-Hamiltonian convention stated in `PROOF.md`.
