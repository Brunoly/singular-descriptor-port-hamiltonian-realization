"""Reproduce exact certificates and adversarial tests. No floating-point proof steps."""
from __future__ import annotations
from pathlib import Path
import json
import random
import sympy as sp
from phdae import (Matrix, kyp, ph_from_Q, split_descriptor, solve_reduced,
                   principal_minors, smt_expr, zero, rectangular_symbols,
                   simplex_determinant_search, is_psd_exact)
from z3_capi import run_smt2

ROOT=Path(__file__).resolve().parent
CERT=ROOT.parent/'checks'
CERT.mkdir(exist_ok=True)


def system(E,A,B,C,D):
    return tuple(Matrix(x) for x in (E,A,B,C,D))


def algebraic_case(n:int,C:Matrix):
    return sp.zeros(n),sp.eye(n),sp.eye(n),C,sp.zeros(n)


def scalar_port(E,A):
    n=E.rows
    return E,A,sp.zeros(n,1),sp.zeros(1,n),sp.zeros(1)


def generic_nra(psd:list[Matrix],equalities:list[sp.Expr],extra:list[str],
                variables:set,timeout:int=20000)->str:
    lines=['(set-logic QF_NRA)',f'(set-option :timeout {timeout})']
    for v in sorted(variables,key=str):lines.append(f'(declare-const {v} Real)')
    for e in equalities:
        if sp.simplify(e)!=0:lines.append(f'(assert (= {smt_expr(sp.expand(e))} 0))')
    seen=set()
    for M in psd:
        for d in principal_minors(M):
            if d==0:continue
            text=smt_expr(d)
            if text not in seen:
                lines.append(f'(assert (>= {text} 0))');seen.add(text)
    lines += ['(assert '+s+')' for s in extra]
    lines.append('(check-sat-using qfnra-nlsat)')
    return '\n'.join(lines)+'\n'


def original_kernel_query(sys,timeout=20000):
    E,A,B,C,D=sys;n=E.rows
    Q=rectangular_symbols('q',n,n)
    M=rectangular_symbols('mult',n,n)
    H=E.T*Q;W=kyp(E,A,B,C,D,Q)
    eq=list(H-H.T)+list(M*Q-A)
    variables=Q.free_symbols|M.free_symbols
    script=generic_nra([H,W],eq,[],variables,timeout)
    return run_smt2(script).strip(),script


def main():
    records=[]
    def record(name,result,detail=''):
        records.append(dict(name=name,result=result,detail=detail))
        print(f'{name}: {result} {detail}',flush=True)
    N2=Matrix([[0,1],[0,0]])
    N3=Matrix([[0,1,0],[0,0,1],[0,0,0]])
    J=Matrix([[0,1],[-1,0]])
    cases={
      'index_one_hidden_damped':(scalar_port(sp.diag(1,0),sp.diag(-1,-1)),'SAT'),
      'index_two_differentiator':((N2,sp.eye(2),Matrix([0,1]),Matrix([[-1,0]]),sp.zeros(1)),'SAT'),
      'hidden_unstable_finite':(scalar_port(sp.diag(1,0),sp.diag(1,1)),'UNSAT'),
      'hidden_index_three':(scalar_port(N3,sp.eye(3)),'UNSAT_INDEX'),
      'forced_singular_algebraic':((sp.zeros(1),sp.eye(1),sp.ones(1),sp.zeros(1),sp.zeros(1)),'UNSAT'),
      'invertible_lossless_algebraic':(algebraic_case(2,J),'SAT'),
      'odd_singular_lossless_algebraic':(algebraic_case(3,sp.diag(J,0)),'UNSAT'),
      'observed_unactuated_no_feedthrough':((sp.diag(1,0),-sp.eye(2),sp.zeros(2,1),Matrix([[1,0]]),sp.zeros(1)),'UNSAT'),
      'allowed_zero_energy_integrator':((sp.diag(1,0),sp.diag(0,-1),Matrix([1,0]),sp.zeros(1,2),sp.zeros(1)),'SAT'),
      'forbidden_zero_jordan_three':(scalar_port(sp.diag(sp.eye(3),0),sp.diag(N3,1)),'UNSAT'),
    }
    E=sp.diag(sp.eye(2),N2,0)
    A=sp.diag(N2,sp.eye(2),1)
    B=Matrix([0,0,0,1,0]);C=Matrix([[0,0,-1,0,0]]);D=sp.zeros(1)
    mixed=(E,A,B,C,D)
    cases['mixed_singular_energy_index_two']=(mixed,'SAT')
    for name,(sys,expected) in cases.items():
        out=solve_reduced(*sys,timeout_ms=25000,output_dir=CERT/name)
        status=out['status']
        assert status==expected,(name,status,expected)
        detail=''
        if status=='SAT':
            detail=f"rank(Q)={out['Q'].rank()}, n={sys[0].rows}"
        record(name,'PASS',status+'; '+detail)

    # Explicit mixed certificate, independent of the solver's selected metric.
    Q=sp.diag(sp.diag(0,1),J,-1)
    ph=ph_from_Q(*mixed,Q)
    assert Q.rank()==4
    record('explicit_mixed_certificate','PASS','all identities and PSD principal minors exact')

    # Strict equivalence and port power preservation under exact rational matrices.
    L=Matrix([[1,1,0,0,0],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1],[0,0,0,0,1]])
    T=Matrix([[1,0,1,0,0],[0,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0],[0,0,0,0,1]])
    U=Matrix([[2]])
    Et,At,Bt,Ct,Dt=L*E*T,L*A*T,L*B*U,U.T*C*T,U.T*D*U
    Qt=L.T.inv()*Q*T
    assert zero(kyp(Et,At,Bt,Ct,Dt,Qt)-sp.diag(T,U).T*kyp(*mixed,Q)*sp.diag(T,U))
    ph_from_Q(Et,At,Bt,Ct,Dt,Qt)
    s=split_descriptor(Et,At,Bt,Ct,Dt)
    assert s.index==2 and s.finite==2 and s.chains_two==1 and s.algebraic_one==1
    record('strict_equivalence_and_rational_split','PASS','nonorthogonal state/equation and scaled power ports')

    # Repeated exact rational changes of coordinates test Fitting decomposition.
    rng=random.Random(20260926)
    for it in range(12):
        Lr=sp.eye(5);Tr=sp.eye(5)
        for _ in range(8):
            i,j=rng.sample(range(5),2);d=rng.choice([-2,-1,1,2])
            Lr[i,:]=Lr[i,:]+d*Lr[j,:]
            i,j=rng.sample(range(5),2);d=rng.choice([-2,-1,1,2])
            Tr[:,i]=Tr[:,i]+d*Tr[:,j]
        ss=split_descriptor(Lr*E*Tr,Lr*A*Tr,Lr*B,C*Tr,D)
        assert (ss.finite,ss.index,ss.chains_two,ss.algebraic_one)==(2,2,1,1)
    record('twelve_exact_scrambled_splits','PASS','all dimensions and identities preserved')

    # Exact rejection of a proposed nonzero finite/infinite energy cross block.
    for nu in (1,2,3):
        N=sp.zeros(nu)
        for i in range(nu-1):N[i,i+1]=1
        Ee=sp.diag(1,N);Aa=sp.diag(Matrix([[-1]]),sp.eye(nu))
        q=rectangular_symbols('q',1+nu,1+nu)
        H=Ee.T*q
        state=-Aa.T*q-q.T*Aa
        cross=sum(q[0,j]**2 for j in range(1,nu+1))
        script=generic_nra([H,state],list(H-H.T),[f'(> {smt_expr(cross)} 0)'],q.free_symbols)
        text=run_smt2(script).strip()
        (CERT/f'cross_block_index_{nu}.smt2').write_text(script)
        (CERT/f'cross_block_index_{nu}.txt').write_text(text+'\n')
        assert text=='unsat',(nu,text)
        record(f'cross_block_counterexample_index_{nu}','PASS','UNSAT for every real Q')

    # Compare the reduced criterion with the unreduced KYP + A=M Q condition.
    comparison_systems=[]
    for i in range(8):
        ff=rng.choice([-1,0,1]);bb=[rng.randint(-1,1) for _ in range(2)]
        cc=[rng.randint(-1,1) for _ in range(2)];dd=rng.choice([0,1])
        comparison_systems.append((sp.diag(1,0),sp.diag(ff,1),Matrix(bb),Matrix([cc]),Matrix([[dd]])))
    for i,sys in enumerate(comparison_systems):
        reduced=solve_reduced(*sys,timeout_ms=15000)['status']
        original,script=original_kernel_query(sys,timeout=15000)
        (CERT/f'original_comparison_{i}.smt2').write_text(script)
        (CERT/f'original_comparison_{i}.txt').write_text(original+'\n')
        assert original in ('sat','unsat') and reduced==original.upper(),(i,reduced,original)
    record('eight_unreduced_vs_reduced_queries','PASS','all exact SAT/UNSAT answers agreed')

    # The polynomial grid can find an invertible convex combination even when
    # every supplied vertex is singular. A separate example is singular everywhere.
    search=simplex_determinant_search([sp.diag(1,0),sp.diag(0,1)])
    assert search is not None and search[0].det()!=0
    assert simplex_determinant_search([Matrix([[0,1,0],[-1,0,0],[0,0,0]]),
                                      Matrix([[0,0,1],[0,0,0],[-1,0,0]])]) is None
    record('determinant_grid','PASS',f'found det={search[0].det()} from singular vertices; odd skew family rejected')

    # A wider behavior-preserving output gauge changes the fixed-coefficient answer.
    E=sp.diag(1,0);A=-sp.eye(2);B=sp.ones(2,1);C=sp.ones(1,2);D=Matrix([[-1]])
    Wg=Matrix([0,-2])
    assert zero(E.T*Wg)
    Cp,Dp=C-Wg.T*A,D-Wg.T*B
    ph_from_Q(E,A,B,Cp,Dp,sp.eye(2))
    assert not is_psd_exact(D+D.T)
    record('output_gauge_scope_example','PASS','strict form obstructed; gauged form certified with Q=I')

    (CERT/'TEST_RESULTS.json').write_text(json.dumps(records,indent=2))
    print(f'\n{len(records)} test groups passed.',flush=True)


if __name__=='__main__':main()
