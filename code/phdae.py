"""Exact, small-dimensional descriptor-to-pH research utilities.

The mathematics is specified in FINAL_PROOF.md.  SymPy is required.  The
optional exact solver uses libz3 through z3_capi.py and deliberately reports
UNKNOWN on solver limits.  The principal-minor encoding is exponential and
is a reference implementation, not a large-scale numerical SDP solver.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
import json
import re
import sympy as sp

Matrix = sp.Matrix


def simplify_matrix(M: Matrix) -> Matrix:
    return M.applyfunc(sp.simplify)


def zero(M: Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in M)


def hcat(cols: list[Matrix], rows: int) -> Matrix:
    return Matrix.hstack(*cols) if cols else sp.zeros(rows, 0)


def principal_minors(M: Matrix) -> list[sp.Expr]:
    if M.rows != M.cols:
        raise ValueError('PSD matrix must be square')
    out: list[sp.Expr] = []
    for k in range(1, M.rows + 1):
        for inds in combinations(range(M.rows), k):
            out.append(sp.factor(M.extract(inds, inds).det()))
    return out


def is_psd_exact(M: Matrix) -> bool:
    if not zero(M - M.T):
        return False
    for minor in principal_minors(M):
        s = sp.simplify(minor)
        if s.is_nonnegative is False:
            return False
        if s.is_nonnegative is not True:
            comparison = sp.ask(sp.Q.nonnegative(s))
            if comparison is not True:
                raise ValueError(f'Exact sign not resolved: {s}')
    return True


def validate_system(E: Matrix, A: Matrix, B: Matrix, C: Matrix, D: Matrix) -> None:
    n, m = E.rows, B.cols
    if E.shape != (n,n) or A.shape != (n,n) or B.rows != n:
        raise ValueError('Invalid state matrix dimensions')
    if C.shape != (m,n) or D.shape != (m,m):
        raise ValueError('Input/output dimensions must be equal')


def kyp(E: Matrix, A: Matrix, B: Matrix, C: Matrix, D: Matrix,
        Q: Matrix) -> Matrix:
    validate_system(E,A,B,C,D)
    L = C.T - Q.T*B
    return Matrix.vstack(Matrix.hstack(-A.T*Q-Q.T*A, L),
                         Matrix.hstack(L.T, D+D.T))


def index_two_basis(N: Matrix) -> Matrix:
    """Rational chain basis [a,b,c] with N b=a and N a=N c=0."""
    h = N.rows
    if not zero(N*N):
        raise ValueError('Index exceeds two')
    aa = N.columnspace()
    bb: list[Matrix] = []
    for a in aa:
        sol, params = N.gauss_jordan_solve(a)
        bb.append(sol.subs({p: 0 for p in params}))
    basis = list(aa)
    cc: list[Matrix] = []
    for c in N.nullspace():
        trial = hcat(basis+[c],h)
        if trial.rank() > len(basis):
            basis.append(c)
            cc.append(c)
    S = hcat(aa+bb+cc,h)
    if S.shape != (h,h) or (h and S.det() == 0):
        raise AssertionError('Failed to construct infinite chain basis')
    return S


@dataclass
class Split:
    E: Matrix
    A: Matrix
    B: Matrix
    C: Matrix
    D: Matrix
    L: Matrix
    T: Matrix
    F: Matrix
    N: Matrix
    finite: int
    chains_two: int
    algebraic_one: int
    index: int
    shift: int


def split_descriptor(E: Matrix,A: Matrix,B: Matrix,C: Matrix,D: Matrix) -> Split:
    """Exact finite/infinite split. Rational input stays rational.

    Resolvent powers are for exact arithmetic ONLY, not floating-point use.
    Reordered generalized Schur/staircase forms are preferable numerically.
    """
    validate_system(E,A,B,C,D)
    n = E.rows
    shifted = None
    alpha = 0
    for alpha in range(n+1):
        cand = alpha*E-A
        if cand.det() != 0:
            shifted = cand
            break
    if shifted is None:
        raise ValueError('Singular pencil: outside theorem hypotheses')
    R = shifted.inv()*E
    Rn = R**n
    uf = Rn.columnspace()
    ui = Rn.nullspace()
    T = hcat(uf+ui,n)
    if n and T.det() == 0:
        raise AssertionError('Fitting decomposition failed')
    r = len(uf); h = n-r
    L = T.inv()*shifted.inv()
    Et = simplify_matrix(L*E*T)
    At = simplify_matrix(L*A*T)
    if not zero(Et[:r,r:]) or not zero(Et[r:,:r]):
        raise AssertionError('Not a deflating decomposition')
    normalizer = sp.diag(Et[:r,:r].inv(), At[r:,r:].inv())
    L = normalizer*L
    Ec = simplify_matrix(L*E*T)
    Ac = simplify_matrix(L*A*T)
    F = Ac[:r,:r]
    N = Ec[r:,r:]
    index = 0 if h == 0 else 1
    if h:
        power = N
        while not zero(power):
            power = power*N
            index += 1
            if index > h:
                raise AssertionError('Infinite block was not nilpotent')
    k = N.rank() if index <= 2 else -1
    a = h-2*k if index <= 2 else -1
    if index <= 2:
        Si = index_two_basis(N)
        S = sp.diag(sp.eye(r),Si)
        T = T*S
        L = S.inv()*L
        Ec = simplify_matrix(L*E*T)
        Ac = simplify_matrix(L*A*T)
        F = Ac[:r,:r]
        N = Ec[r:,r:]
    Bc, Cc = simplify_matrix(L*B), simplify_matrix(C*T)
    assert zero(Ec-sp.diag(sp.eye(r),N))
    assert zero(Ac-sp.diag(F,sp.eye(h)))
    return Split(Ec,Ac,Bc,Cc,D,L,T,F,N,r,k,a,index,alpha)


def ph_from_Q(E: Matrix,A: Matrix,B: Matrix,C: Matrix,D: Matrix,Q: Matrix,
              strong_dissipation: bool = True) -> dict[str,Matrix]:
    """Factor a certified KYP variable without inverting its singular part."""
    W = kyp(E,A,B,C,D,Q)
    if not is_psd_exact(E.T*Q) or not is_psd_exact(W):
        raise ValueError('Not a KYP certificate')
    Qp = Q.pinv()
    if not zero(A*Qp*Q-A) or not zero(C*Qp*Q-C):
        raise ValueError('Kernel compatibility fails')
    n,m = E.rows,B.cols
    if strong_dissipation:
        QQ = sp.diag(Q,sp.eye(m))
        AA = Matrix.vstack(Matrix.hstack(A,B), Matrix.hstack(-C,-D))
        theta0 = AA*QQ.pinv()
        Pi = QQ*QQ.pinv()
        theta = simplify_matrix(theta0-Pi*theta0.T*(sp.eye(n+m)-Pi))
        skew = (theta-theta.T)/2
        diss = -(theta+theta.T)/2
        J,R = skew[:n,:n],diss[:n,:n]
        G,P = skew[:n,n:],diss[:n,n:]
        S,Nport = diss[n:,n:],-skew[n:,n:]
    else:
        M, Fo = A*Qp,C*Qp
        J,R = (M-M.T)/2,-(M+M.T)/2
        G,P = (B+Fo.T)/2,(Fo.T-B)/2
        S,Nport = (D+D.T)/2,(D-D.T)/2
    ans = {name:simplify_matrix(M) for name,M in
           dict(E=E,Q=Q,J=J,R=R,G=G,P=P,S=S,N=Nport).items()}
    assert zero((J-R)*Q-A)
    assert zero(G-P-B)
    assert zero((G+P).T*Q-C)
    assert zero(S+Nport-D)
    Dw = Matrix.vstack(Matrix.hstack(Q.T*R*Q,Q.T*P),
                        Matrix.hstack(P.T*Q,S))
    assert zero(2*Dw-W)
    assert is_psd_exact(Dw)
    if strong_dissipation:
        diss = Matrix.vstack(Matrix.hstack(R,P),Matrix.hstack(P.T,S))
        assert is_psd_exact(diss)
    return ans


def symmetric_symbols(prefix: str,n: int) -> Matrix:
    syms = {(i,j):sp.Symbol(f'{prefix}_{i}_{j}',real=True)
            for i in range(n) for j in range(i,n)}
    return Matrix(n,n,lambda i,j:syms[min(i,j),max(i,j)])


def rectangular_symbols(prefix: str,r: int,c: int) -> Matrix:
    return Matrix(r,c,lambda i,j:sp.Symbol(f'{prefix}_{i}_{j}',real=True))


@dataclass
class ReducedProblem:
    split: Split
    Q: Matrix
    X: Matrix
    K: Matrix
    Z00: Matrix
    epsilon: sp.Symbol
    psd: list[Matrix]
    variables: list[sp.Symbol]


def reduced_problem(s: Split) -> ReducedProblem:
    if s.index > 2:
        raise ValueError('Index obstruction')
    r,k,a = s.finite,s.chains_two,s.algebraic_one
    h=2*k+a
    X=symmetric_symbols('X',r)
    K=symmetric_symbols('K',k)
    Y2=rectangular_symbols('Y2',k,r)
    Y0=rectangular_symbols('Y0',a,r)
    Z22=rectangular_symbols('Z22',k,k)
    Z20=rectangular_symbols('Z20',k,a)
    Z02=rectangular_symbols('Z02',a,k)
    Z00=rectangular_symbols('Z00',a,a)
    Y=Matrix.vstack(sp.zeros(k,r),Y2,Y0)
    Z=Matrix.vstack(Matrix.hstack(sp.zeros(k,k),K,sp.zeros(k,a)),
                    Matrix.hstack(-K,Z22,Z20),
                    Matrix.hstack(sp.zeros(a,k),Z02,Z00))
    Q=Matrix.vstack(Matrix.hstack(X,sp.zeros(r,h)),Matrix.hstack(Y,Z))
    eps=sp.Symbol('epsilon',real=True)
    psd=[X-eps*s.F.T*s.F,K-eps*sp.eye(k),
         kyp(s.E,s.A,s.B,s.C,s.D,Q)]
    syms=set([eps])
    for M in psd:
        syms.update(M.free_symbols)
    syms.update(Q.free_symbols)
    return ReducedProblem(s,Q,X,K,Z00,eps,psd,sorted(syms,key=str))


def smt_expr(e: sp.Expr) -> str:
    e=sp.sympify(e)
    if e.is_Symbol:
        return str(e)
    if e.is_Rational:
        p,q=int(e.p),int(e.q)
        ps=str(abs(p)) if p>=0 else f'(- {abs(p)})'
        return ps if q==1 else f'(/ {ps} {q})'
    if e.is_Add:
        return '(+ '+' '.join(smt_expr(x) for x in e.args)+')'
    if e.is_Mul:
        return '(* '+' '.join(smt_expr(x) for x in e.args)+')'
    if e.is_Pow and e.exp.is_Integer and e.exp>=0:
        n=int(e.exp)
        if n==0:return '1'
        return '(* '+' '.join([smt_expr(e.base)]*n)+')' if n>1 else smt_expr(e.base)
    raise ValueError(f'Non-polynomial or non-rational SMT coefficient: {e}')


def export_smt2(p: ReducedProblem, timeout_ms: int|None=None,
                extra: list[str]|None=None, determinant: bool=True) -> str:
    lines=['(set-logic QF_NRA)','(set-option :produce-models true)']
    if timeout_ms is not None:
        lines.append(f'(set-option :timeout {int(timeout_ms)})')
    for v in p.variables:
        lines.append(f'(declare-const {v} Real)')
    lines.append(f'(assert (> {p.epsilon} 0))')
    seen=set()
    for M in p.psd:
        assert zero(M-M.T)
        for d in principal_minors(M):
            if d==0:continue
            s=smt_expr(d)
            if s not in seen:
                seen.add(s)
                lines.append(f'(assert (>= {s} 0))')
    if determinant and p.Z00.rows:
        lines.append(f'(assert (not (= {smt_expr(sp.expand(p.Z00.det()))} 0)))')
    for s in extra or []:
        lines.append(f'(assert {s})')
    lines += ['(check-sat-using qfnra-nlsat)',
              '(get-value ('+' '.join(map(str,p.variables))+'))']
    return '\n'.join(lines)+'\n'


def parse_sexps(text: str) -> list:
    tokens=re.findall(r'\(|\)|"(?:[^"\\]|\\.)*"|[^\s()]+',text)
    i=0
    def parse():
        nonlocal i
        if i>=len(tokens):raise ValueError('Unexpected end of SMT output')
        tok=tokens[i];i+=1
        if tok!='(':return tok
        seq=[]
        while i<len(tokens) and tokens[i]!=')':seq.append(parse())
        if i>=len(tokens):raise ValueError('Unbalanced SMT output')
        i+=1
        return seq
    out=[]
    while i<len(tokens):out.append(parse())
    return out


def exact_smt_value(v, env: dict|None=None) -> sp.Expr:
    env=env or {}
    if isinstance(v,str):
        return env[v] if v in env else sp.Rational(v)
    op,*args=v
    if op=='root-obj':
        t=sp.Symbol('x')
        poly=exact_smt_value(args[0],{'x':t})
        roots=sp.Poly(poly,t).real_roots(radicals=False)
        return roots[int(args[1])-1]
    vals=[exact_smt_value(x,env) for x in args]
    if op=='+':return sp.Add(*vals)
    if op=='*':return sp.Mul(*vals)
    if op=='-':return -vals[0] if len(vals)==1 else vals[0]-vals[1]
    if op=='/':return vals[0]/vals[1]
    if op=='^':return vals[0]**vals[1]
    raise ValueError(f'Unsupported exact SMT value: {v}')


def solve_reduced(E: Matrix,A: Matrix,B: Matrix,C: Matrix,D: Matrix,
                  timeout_ms: int|None=20000, output_dir: str|Path|None=None) -> dict:
    """Reference exact solver for rational input, with explicit UNKNOWN status.

    With timeout_ms=None no solver time limit is imposed.  Correctness of the
    returned positive result is rechecked independently by exact substitution.
    An UNSAT answer relies on the SMT engine; the main theorem does not.
    """
    from z3_capi import run_smt2
    if output_dir:
        path=Path(output_dir);path.mkdir(parents=True,exist_ok=True)
        original={name:[[str(x) for x in M.row(i)] for i in range(M.rows)]
                  for name,M in dict(E=E,A=A,B=B,C=C,D=D).items()}
        (path/'system.json').write_text(json.dumps(original,indent=2)+'\n')
    s=split_descriptor(E,A,B,C,D)
    if s.index>2:
        if output_dir:
            (Path(output_dir)/'index_obstruction.json').write_text(
                json.dumps({'infinite_index':s.index,'conclusion':'UNSAT_INDEX'},indent=2)+'\n')
        return dict(status='UNSAT_INDEX',index=s.index,split=s)
    p=reduced_problem(s)
    script=export_smt2(p,timeout_ms=timeout_ms)
    out=run_smt2(script)
    if output_dir:
        path=Path(output_dir);path.mkdir(parents=True,exist_ok=True)
        (path/'problem.smt2').write_text(script)
        (path/'solver_output.txt').write_text(out)
    sexps=parse_sexps(out)
    if not sexps or sexps[0] not in ('sat','unsat','unknown'):
        raise RuntimeError(f'Invalid solver response: {out}')
    status=sexps[0]
    if status!='sat':return dict(status=status.upper(),raw=out,split=s)
    values={str(v):v for v in p.variables}
    if len(sexps)<2:raise RuntimeError('Missing model')
    subs={values[name]:exact_smt_value(val) for name,val in sexps[1]}
    Qc=simplify_matrix(p.Q.subs(subs))
    for M in p.psd:
        if not is_psd_exact(simplify_matrix(M.subs(subs))):
            raise AssertionError('Model failed exact PSD check')
    eps=sp.simplify(subs[p.epsilon])
    if eps.is_positive is not True:raise AssertionError('Nonpositive margin')
    Z00=simplify_matrix(p.Z00.subs(subs))
    if Z00.det()==0:raise AssertionError('Singular algebraic certificate')
    Q=simplify_matrix(s.L.T*Qc*s.T.inv())
    factors=ph_from_Q(E,A,B,C,D,Q)
    ans=dict(status='SAT',Q=Q,Q_canonical=Qc,epsilon=eps,
             factors=factors,split=s,raw=out)
    if output_dir:
        serial={name:[[str(x) for x in M.row(i)] for i in range(M.rows)]
                for name,M in factors.items()}
        serial['epsilon']=str(eps)
        (Path(output_dir)/'certificate.json').write_text(json.dumps(serial,indent=2))
    return ans


def simplex_determinant_search(matrices: list[Matrix]) -> tuple[Matrix,list[sp.Rational]]|None:
    """Finite interpolation grid. Assumes supplied points span the feasible affine hull.

    This function does NOT certify that spanning assumption: the affine-hull
    oracle procedure is given separately in FINAL_PROOF.md.
    """
    if not matrices:raise ValueError('At least one feasible point is required')
    a=matrices[0].rows
    if any(M.shape!=(a,a) for M in matrices):raise ValueError('Inconsistent shapes')
    # Keep an affine basis, recording indices to report weights on original points.
    inds=[0]; cols=[]
    for j,M in enumerate(matrices[1:],1):
        v=Matrix(list(M-matrices[0]))
        if hcat(cols+[v],a*a).rank()>len(cols):
            cols.append(v);inds.append(j)
    d=len(inds)-1
    denom=d*(a+1)+1
    grid=product(range(1,a+2),repeat=d)
    for numbers in grid:
        tail=[sp.Rational(j,denom) for j in numbers]
        lam=[1-sum(tail)]+tail
        candidate=sp.zeros(a,a)
        for w,i in zip(lam,inds):candidate+=w*matrices[i]
        if candidate.det()!=0:
            weights=[sp.Rational(0)]*len(matrices)
            for w,i in zip(lam,inds):weights[i]=w
            return candidate,weights
    return None
