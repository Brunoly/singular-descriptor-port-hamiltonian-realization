"""Independent checker for the bundled rational pH certificates.

Requires only SymPy, not libz3 or phdae.py. Inputs are parsed strictly as
rational numbers, not executable expressions. General algebraic solver models
are checked within phdae.py; this separate checker targets the saved rational
regression certificates.
"""
from __future__ import annotations
from itertools import combinations
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parent


def matrix(data: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry) for entry in row] for row in data])


def require_zero(M: sp.Matrix, label: str) -> None:
    if any(entry != 0 for entry in M):
        raise AssertionError(f'Nonzero identity residual: {label}')


def require_psd(M: sp.Matrix, label: str) -> None:
    require_zero(M-M.T, label+' symmetry')
    for size in range(1,M.rows+1):
        for inds in combinations(range(M.rows),size):
            value=M.extract(inds,inds).det()
            if value < 0:
                raise AssertionError(f'Negative principal minor: {label}, {inds}: {value}')


def check(path: Path) -> None:
    raw=json.loads(path.read_text())
    system=json.loads((path.parent/'system.json').read_text())
    E,A,B,C,D=(matrix(system[name]) for name in ('E','A','B','C','D'))
    Q,J,R,G,P,S,N=(matrix(raw[name]) for name in ('Q','J','R','G','P','S','N'))
    require_zero(E-matrix(raw['E']),'descriptor coefficient')
    require_zero(J+J.T,'J skew')
    require_zero(N+N.T,'N skew')
    require_zero(R-R.T,'R symmetric')
    require_zero(S-S.T,'S symmetric')
    require_zero((J-R)*Q-A,'A factorization')
    require_zero(G-P-B,'B factorization')
    require_zero((G+P).T*Q-C,'C factorization')
    require_zero(S+N-D,'D factorization')
    require_psd(E.T*Q,'energy')
    diss=sp.Matrix.vstack(sp.Matrix.hstack(R,P),sp.Matrix.hstack(P.T,S))
    weighted=sp.diag(Q,sp.eye(B.cols)).T*diss*sp.diag(Q,sp.eye(B.cols))
    require_psd(weighted,'weighted dissipation')
    require_psd(diss,'strong unweighted dissipation')
    print(path.parent.name+': independent rational certificate PASS')


def main() -> None:
    paths=sorted((ROOT.parent/'checks').glob('*/certificate.json'))
    if not paths:
        raise RuntimeError('No saved certificates found')
    for path in paths:
        check(path)
    print(f'\n{len(paths)} independent rational certificate checks passed.')


if __name__=='__main__':
    main()
