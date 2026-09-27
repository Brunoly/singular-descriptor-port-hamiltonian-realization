# Exact realization of regular descriptor systems in port-Hamiltonian form

**Date:** September 26, 2026.  
**Status:** Complete mathematical characterization and finite exact construction for the explicitly defined equivalence classes below. Research draft; not independently peer-reviewed or formally verified. Historical priority is not established by the literature search.

## Abstract

For an arbitrary real regular descriptor realization, including nonminimal realizations and singular symmetric feedthrough, we reduce port-Hamiltonian realizability to a convex system of linear matrix inequalities and one determinant condition confined to the index-one algebraic block. The determinant condition is decided by finitely many convex feasibility problems followed by a finite polynomial interpolation test. No controllability, observability, finite stability, index-one, or nonsingular-energy assumption is imposed in advance.

The essential structural lemma is that every global descriptor KYP solution in finite/infinite deflating coordinates is block lower triangular. Its finite block is positive semidefinite, its infinite block satisfies a nilpotent dissipation identity, and its off-diagonal block annihilates the finite energy kernel. Consequently the previously known kernel compatibility condition becomes an ordinary LMI on the finite block plus invertibility of the infinite block. This also gives an independent proof that a regular realization of the specified pH form cannot have index greater than two. For index two, the remaining nonsingularity question lies entirely in the size-one infinite blocks.

The algorithm is finite in exact arithmetic for rational or real-algebraic input. It is not claimed to be polynomial-time or uniformly numerically stable. An exact small-dimensional reference implementation, scripts, solver queries, outputs, and positive certificates accompany the proof.

---

## 1. Specification, notation, and scope

Consider

$$
 E\dot x=Ax+Bu,\qquad y=Cx+Du,
 \qquad\text{(1.1)}
$$

where $E,A\in\mathbb R^{n\times n}$, $B\in\mathbb R^{n\times m}$, $C\in\mathbb R^{m\times n}$, and $D\in\mathbb R^{m\times m}$. Assume only

$$
 \det(sE-A)\not\equiv0.
 \qquad\text{(1.2)}
$$

The motivating problem assumes singular $E$, but the theorem also includes invertible $E$. All statements are over the reals; no reduction to complex coordinates is necessary. Transpose is denoted by $T$, the Moore–Penrose inverse by $\dagger$, and positive semidefiniteness by $\succeq0$.

A **pH representation** means the same coefficient realization, or a realization equivalent in the specified sense, written as

$$
 E\dot x=(J-R)Qx+(G-P)u,
 \qquad y=(G+P)^TQx+(S+N_p)u,
 \qquad\text{(1.3)}
$$

with

$$
 J=-J^T,\quad R=R^T,\quad S=S^T,\quad N_p=-N_p^T,
 \qquad\text{(1.4)}
$$

$$
 E^TQ=Q^TE\succeq0,
 \qquad
 \begin{bmatrix}Q^TRQ&Q^TP\\P^TQ&S\end{bmatrix}\succeq0.
 \qquad\text{(1.5)}
$$

The matrix $Q$ may be singular. The subscript in $N_p$ distinguishes the skew port matrix from the nilpotent infinite-eigenvalue matrix $N$.

The default equivalence class consists of **constant invertible state transformations, invertible equation transformations, and power-preserving invertible port transformations**, as defined in Section 3. It does not delete nonminimal states, enlarge the state space, apply feedback, or replace a system merely by another one with the same transfer function. Section 12 gives a second complete theorem for a larger, explicitly defined class allowing algebraic output gauges.

For the trajectory interpretation, take differentiable state trajectories and inputs/outputs for which (1.1) is meaningful. Only admissible trajectories are included: no assumption is made that arbitrary nonsmooth inputs admit classical descriptor solutions. The transformations give bijections of these trajectories and consistency conditions. Constant invertible transformations also respect distributional versions of the state equations; the pointwise supply product is used here only for classical trajectories.

Zero-dimensional subblocks are allowed throughout. A zero-by-zero identity or PSD constraint is vacuous, and $\det(0_{0\times0})=1$.

### What is, and is not, claimed

The result is a full iff theorem plus a finite construction for the above target and equivalences. It is not a claim that all passive or positive-real realizations are strictly equivalent to pH realizations. It is not a first-ever descriptor pH realization theorem. The basic KYP-plus-kernel theorem is already known [R2], and general positive-real transfer-function realization procedures are already known [R3]. The contribution developed in this document is the explicit reduction and finite determinant decision described below. An exhaustive priority determination would require additional independent bibliographic and expert review.

## 2. Literature boundary checked

The foundational descriptor pH framework is [R1]. Cherifi–Gernandt–Hinsen [R2, Proposition 11] already establish the equivalence between pH representation and a global KYP solution satisfying $\ker Q\subseteq\ker A\cap\ker C$. Their Remark 12 discusses a dissipative completion for singular $Q$. We use and independently prove the relevant factorization below; neither fact is presented as new.

The 2024 Chu–Mehrmann preprint arXiv:2408.14115 has a published version [R3], *Automatica* 180 (2025), 112456. Its Theorem 7 gives a fixed-coefficient characterization under complete controllability and complete observability: for positive-real systems, pH representation is equivalent to $D+D^T\succeq0$. Corollary 8 supplies the corresponding KYP equivalence. Theorem 10 uses an algebraic output correction. Corollary 11 and the subsequent algorithm handle general positive-real transfer-function realization after controllability/observability reduction. That reduction need not retain all states of the supplied realization.

Beattie–Mehrmann–Xu [R4, Remark 30] explicitly identify the extension of their nonminimal standard-system constructive framework to singular $E$ as open. Remark 31 notes their restriction to positive-definite energy metrics. Merely repeating [R2] would therefore not address that constructive extension.

The index-at-most-two restriction for regular dissipative-Hamiltonian pencils is established in earlier literature [R5]; Section 5 proves the needed restriction directly for the exact convention (1.3)–(1.5). It is not a blanket assertion about every broader constrained, rectangular, or nonregular pH formalism.

The June 26, 2026 regularization paper [R6] concerns proportional/derivative output feedback of systems already in pH form. It does not replace the state-preserving realization question here. The audit checked the named primary papers and related current results; it does not prove the absence of every unindexed or differently formulated equivalent result.

## 3. Exact equivalence and invariance

Let $L,T\in GL_n(\mathbb R)$ and $V\in GL_m(\mathbb R)$. Set

$$
 x=Tz,\qquad u=Vv,\qquad \eta=V^Ty,
$$

and left-multiply the state equations by $L$. The transformed realization is

$$
 \widehat E=LET,\quad \widehat A=LAT,\quad
 \widehat B=LBV,\quad \widehat C=V^TCT,\quad
 \widehat D=V^TDV.
 \qquad\text{(3.1)}
$$

All three matrices are invertible, so this is a trajectory bijection without elimination or additional solutions. Moreover

$$
 \eta^Tv=y^TVv=y^Tu.
 \qquad\text{(3.2)}
$$

Define the global KYP matrix

$$
 \mathcal W_\Sigma(Q)=
 \begin{bmatrix}
 -A^TQ-Q^TA&C^T-Q^TB\\
 C-B^TQ&D+D^T
 \end{bmatrix}.
 \qquad\text{(3.3)}
$$

Under (3.1), take

$$
 \widehat Q=L^{-T}QT.
 \qquad\text{(3.4)}
$$

Direct multiplication gives

$$
 \widehat E^T\widehat Q=T^TE^TQT,
 \qquad
 \mathcal W_{\widehat\Sigma}(\widehat Q)
 =\mathrm{diag}(T,V)^T\mathcal W_\Sigma(Q)
   \mathrm{diag}(T,V).
 \qquad\text{(3.5)}
$$

Also $\ker\widehat Q=T^{-1}\ker Q$, and therefore

$$
 \ker Q\subseteq\ker A\cap\ker C
 \iff
 \ker\widehat Q\subseteq\ker\widehat A\cap\ker\widehat C.
 \qquad\text{(3.6)}
$$

For an existing pH factorization, the transformed factors are

$$
 \widehat J=LJL^T,\quad\widehat R=LRL^T,\quad
 \widehat G=LGV,\quad\widehat P=LPV,
$$

$$
 \widehat S=V^TSV,\qquad\widehat N_p=V^TN_pV.
 \qquad\text{(3.7)}
$$

These identities prove invariance of pH realizability under precisely the transformations claimed. In particular, finding a pH realization anywhere in this equivalence class is equivalent to factoring the original realization itself.

## 4. Elementary lemmas and the known kernel criterion

### Lemma 4.1: a zero quadratic form in a PSD matrix

If $M=M^T\succeq0$ and $v^TMv=0$, then $Mv=0$.

**Proof.** Diagonalize $M$ orthogonally, or write $M=H^TH$. Then $v^TMv=\|Hv\|^2=0$, so $Mv=H^THv=0$. ∎

### Lemma 4.2: finite-dimensional domination

For $X=X^T\succeq0$ and an arbitrary matrix $F$ with the same number of columns,

$$
 \ker X\subseteq\ker F
 \iff \exists\epsilon>0:\quad X\succeq\epsilon F^TF.
 \qquad\text{(4.1)}
$$

**Proof.** The reverse implication follows by testing the inequality on $\ker X$. For the forward implication, decompose $v=v_0+v_1$ orthogonally with $v_0\in\ker X$ and $v_1\perp\ker X$. If the complement is nonzero, let $\lambda>0$ be the smallest eigenvalue of the restriction of $X$ to that complement. Then

$$
 v^TXv\ge\lambda\|v_1\|^2,\qquad
 \|Fv\|^2=\|Fv_1\|^2\le\|F\|^2\|v_1\|^2.
$$

When $F\ne0$, choose $0<\epsilon\le\lambda/\|F\|^2$. If $F=0$, any positive $\epsilon$ works. If the complement is zero, the assumed kernel inclusion already forces $F=0$. ∎

### Lemma 4.3: KYP automatically annihilates output on its kernel

If $\mathcal W_\Sigma(Q)\succeq0$, then

$$
 \ker Q\subseteq\ker C,
 \qquad Q^TA\ker Q=\{0\}.
 \qquad\text{(4.2)}
$$

**Proof.** For $Qv=0$, the quadratic form of $\mathcal W_\Sigma(Q)$ at $(v,0)$ is zero. Lemma 4.1 makes the entire product zero. Its state part is $-Q^TAv$, and its port part is $Cv$. ∎

### Proposition 4.4: KYP with factorization compatibility

The realization has a pH factorization (1.3)–(1.5) if and only if there is a real $Q$ such that

$$
 E^TQ=Q^TE\succeq0,\qquad
 \mathcal W_\Sigma(Q)\succeq0,\qquad
 \ker Q\subseteq\ker A.
 \qquad\text{(4.3)}
$$

This is the known kernel criterion [R2, Proposition 11], using Lemma 4.3 to omit the redundant output-kernel inclusion.

**Proof.** For a pH factorization, $A=(J-R)Q$, so $A\ker Q=0$. Substitution gives

$$
 \mathcal W_\Sigma(Q)=2
 \begin{bmatrix}Q^TRQ&Q^TP\\P^TQ&S\end{bmatrix}\succeq0.
 \qquad\text{(4.4)}
$$

Conversely, (4.3) and Lemma 4.3 give $A=AQ^\dagger Q$ and $C=CQ^\dagger Q$. Set

$$
 M=AQ^\dagger,\qquad F_o=CQ^\dagger,
$$

$$
 J=\tfrac12(M-M^T),\quad R=-\tfrac12(M+M^T),
$$

$$
 G=\tfrac12(B+F_o^T),\quad P=\tfrac12(F_o^T-B),
$$

$$
 S=\tfrac12(D+D^T),\quad N_p=\tfrac12(D-D^T).
 \qquad\text{(4.5)}
$$

The factorization identities follow immediately. Its weighted dissipation matrix equals $\mathcal W_\Sigma(Q)/2$, proving (1.5). No inverse of a singular $Q$ is taken. ∎

The storage $\mathcal H(x)=\frac12x^TE^TQx$ then satisfies, along every classical trajectory,

$$
 y^Tu-\dot{\mathcal H}
 =\tfrac12\begin{bmatrix}x\\u\end{bmatrix}^{T}
   \mathcal W_\Sigma(Q)
   \begin{bmatrix}x\\u\end{bmatrix}\ge0.
 \qquad\text{(4.6)}
$$

No converse from passivity or positive realness is assumed.

## 5. The finite/infinite separation lemma

Use real strict equivalence to obtain

$$
 E_c=\begin{bmatrix}I_f&0\\0&N\end{bmatrix},\qquad
 A_c=\begin{bmatrix}F&0\\0&I_h\end{bmatrix},
 \qquad\text{(5.1)}
$$

where $N$ is nilpotent, $f+h=n$, and $F$ is an arbitrary real matrix. Section 10 proves a finite exact construction of these coordinates without requiring Jordan form for $F$. Partition $B_c=(B_f^T,B_\infty^T)^T$ and $C_c=(C_f,C_\infty)$.

### Lemma 5.1: nilpotent Lyapunov identity

If $N$ is nilpotent, $H=H^T\succeq0$, and

$$
 N^TH+HN\preceq0,
$$

then $HN=0$.

**Proof.** For every real vector $v$,

$$
 g_v(t)=\|H^{1/2}e^{tN}v\|^2
$$

is nonnegative and has derivative

$$
 g_v'(t)=v^Te^{tN^T}(N^TH+HN)e^{tN}v\le0
$$

for $t\ge0$. Since $N$ is nilpotent, the vector $H^{1/2}e^{tN}v$ is a polynomial in $t$. Its norm is bounded for $t\ge0$, so every component polynomial is constant. Differentiating at zero gives $H^{1/2}Nv=0$. This holds for every $v$, and hence $HN=0$. ∎

### Theorem 5.2: structure of every global KYP solution

Every $Q$ satisfying

$$
 E_c^TQ=Q^TE_c\succeq0,\qquad \mathcal W_c(Q)\succeq0
$$

has the form

$$
 Q=\begin{bmatrix}X&0\\Y&Z\end{bmatrix},
 \qquad\text{(5.2)}
$$

where

$$
 X=X^T\succeq0,\quad N^TY=0,\quad
 H=N^TZ=Z^TN\succeq0,
 \qquad\text{(5.3)}
$$

$$
 ZN=-N^TZ=-H,\qquad HN=0,\qquad ZN^2=0.
 \qquad\text{(5.4)}
$$

Furthermore,

$$
 Y\ker X=0,\qquad XF\ker X=0,\qquad C_f\ker X=0.
 \qquad\text{(5.5)}
$$

**Proof.** Initially write

$$
 Q=\begin{bmatrix}X&U\\Y&Z\end{bmatrix}.
$$

Energy symmetry and principal-block positivity give

$$
 X=X^T\succeq0,\quad U=Y^TN,\quad
 H=N^TZ=Z^TN\succeq0.
 \qquad\text{(5.6)}
$$

The infinite principal block of $\mathcal W_c(Q)$ is $-Z-Z^T$, so $Z+Z^T\preceq0$. Consequently

$$
 N^TH+HN=N^T(Z+Z^T)N\preceq0.
$$

Lemma 5.1 gives $HN=0$. Thus the quadratic form of $-Z-Z^T$ vanishes on $\mathrm{im}N$. Lemma 4.1 yields

$$
 (Z+Z^T)N=0.
$$

It follows that $ZN=-Z^TN=-H$, and then $ZN^2=0$.

Now apply $\mathcal W_c(Q)$ to a vector whose finite, infinite, and port components are $(0,Nw,0)$. Its quadratic form is zero, so its entire image is zero. The finite component gives

$$
 (F^TU+Y^T)N=0.
$$

Since $Y^TN=U$, this is

$$
 U=-F^TUN.
$$

Iteration gives $U=(-F^T)^jUN^j$. Choosing $j$ with $N^j=0$ proves $U=0$, and (5.6) then gives $N^TY=0$.

Finally take $v\in\ker X$ and apply the PSD zero-vector argument to $(v,0,0)$. Its finite component is $-XFv$, its infinite component is $-Yv$, and, after these vanish, its port component is $C_fv$. This proves (5.5). ∎

### Corollary 5.3: elimination of the kernel condition

For a KYP solution in (5.2),

$$
 \ker Q\subseteq\ker A_c
 \iff Z\text{ is invertible and }\ker X\subseteq\ker F.
 \qquad\text{(5.7)}
$$

**Proof.** If $Zw=0$, then $Q(0,w)=0$, whereas $A_c(0,w)=(0,w)$. Kernel compatibility therefore forces $w=0$, so $Z$ is invertible. For $v\in\ker X$, (5.5) implies $Yv=0$, hence $(v,0)\in\ker Q$ and $Fv=0$.

Conversely, if $Z$ is invertible, then (5.5) shows

$$
 \ker Q=\ker X\oplus\{0\}.
$$

The claimed compatibility follows from $F\ker X=0$. ∎

Combining Lemma 4.2 and Corollary 5.3 gives the first reduced iff:

$$
 \boxed{\quad
 \text{pH realizability}
 \iff
 \exists Q,\epsilon>0:
 \begin{cases}
 Q\text{ satisfies the global KYP conditions},\\
 Q=\begin{bmatrix}X&0\\Y&Z\end{bmatrix},\\
 X\succeq\epsilon F^TF,\\
 \det Z\ne0.
 \end{cases}\quad}
 \qquad\text{(5.8)}
$$

Because $ZN^2=0$, an invertible $Z$ forces $N^2=0$. Thus:

> **Every regular descriptor realization of the specified pH form has infinite index at most two. A regular supplied realization with index greater than two is not strictly equivalent to that form, even when all its high-index states are invisible in the transfer function.**

This rejection is an exact obstruction, not an unhandled case.

## 6. Main theorem: the only residual determinant is algebraic index one

Assume $N^2=0$, and put

$$
 k=\mathrm{rank}N,\qquad a=h-2k.
$$

Choose chain coordinates of sizes $(k,k,a)$ with

$$
 N=\begin{bmatrix}0&I_k&0\\0&0&0\\0&0&0\end{bmatrix}.
 \qquad\text{(6.1)}
$$

There are $k$ infinite Jordan blocks of size two and $a$ of size one. No Jordan form for finite dynamics is used.

Introduce real variables

$$
 X=X^T\in\mathbb R^{f\times f},\quad
 K=K^T\in\mathbb R^{k\times k},
$$

$$
 Y_2\in\mathbb R^{k\times f},\quad Y_0\in\mathbb R^{a\times f},
$$

and arbitrary blocks $Z_{22},Z_{20},Z_{02},Z_{00}$ of sizes $k\times k,k\times a,a\times k,a\times a$. Form

$$
 Q_c=
 \begin{bmatrix}
 X&0&0&0\\
 0&0&K&0\\
 Y_2&-K&Z_{22}&Z_{20}\\
 Y_0&0&Z_{02}&Z_{00}
 \end{bmatrix}.
 \qquad\text{(6.2)}
$$

Let $\mathscr F$ be the convex set of these blocks and a scalar $\epsilon$ satisfying

$$
 \epsilon>0,\qquad
 X-\epsilon F^TF\succeq0,\qquad
 K-\epsilon I_k\succeq0,\qquad
 \mathcal W_c(Q_c)\succeq0.
 \qquad\text{(6.3)}
$$

Every matrix inequality in (6.3) is affine in the unknowns.

### Theorem 6.1: full exact characterization

For any real regular descriptor realization (1.1), the following statements are equivalent:

**(i)** It is equivalent under (3.1) to a pH realization satisfying (1.3)–(1.5).

**(ii)** Its infinite nilpotent block satisfies $N^2=0$, and the set $\mathscr F$ contains a tuple with

$$
 \det Z_{00}\ne0.
 \qquad\text{(6.4)}
$$

**(iii)** Its infinite index is at most two, $\mathscr F\ne\varnothing$, and the determinant polynomial is not identically zero on the affine hull of the projection

$$
 \mathscr Z=\{Z_{00}:\text{some feasible tuple in }\mathscr F\}.
 \qquad\text{(6.5)}
$$

The equivalence includes arbitrary uncontrollable and unobservable finite and algebraic states, arbitrary $D+D^T$, and singular $Q$. A successful tuple constructs pH factors of the original realization by (3.4) inverted and (4.5). The exact finite algorithm of Section 8 decides (ii)–(iii) and produces such a tuple when one exists.

### Proof of the block structure and determinant identity

Write the blocks of $Z$ according to (6.1). The matrix $N^TZ$ consists of the first block row of $Z$ placed in its second block row, with its other block rows zero. Its symmetry and PSD property force

$$
 Z_{11}=0,\quad Z_{13}=0,\quad Z_{12}=K=K^T\succeq0.
$$

The identity $ZN=-N^TZ$ in Theorem 5.2 then forces $Z_{21}=-K$ and $Z_{31}=0$. Also $N^TY=0$ forces the first block row of $Y$ to be zero. Therefore every KYP solution has precisely the shape (6.2), with

$$
 E_c^TQ_c=\mathrm{diag}(X,0_k,K,0_a).
 \qquad\text{(6.6)}
$$

In these coordinates,

$$
 Z=
 \begin{bmatrix}0&K&0\\-K&Z_{22}&Z_{20}\\0&Z_{02}&Z_{00}\end{bmatrix},
 \qquad
 \det Z=(\det K)^2\det Z_{00}.
 \qquad\text{(6.7)}
$$

For completeness, when $K$ is invertible the upper-left two-by-two block has determinant $(\det K)^2$. Solving its linear system on the right-hand side $(0,Z_{20})^T$ gives $(-K^{-1}Z_{20},0)^T$, so the Schur complement is exactly $Z_{00}$. This proves (6.7) on the dense set of invertible $K$; both sides are polynomial in all entries, so it holds identically, including singular $K$.

As $K\succeq0$, invertibility of $Z$ is equivalent to $K\succ0$ and $\det Z_{00}\ne0$. Corollary 5.3 and Lemma 4.2 therefore give exactly (6.3)–(6.4). A common positive $\epsilon$ can be chosen by reducing the separate positive bounds for $X$ and $K$. Conversely, (6.3) gives $X,K\succeq0$, so (6.6) supplies the energy condition; (6.7) and (6.4) give invertible $Z$; then Corollary 5.3 and Proposition 4.4 construct pH factors. This proves (i)⇔(ii).

Section 8 proves (ii)⇔(iii), including a finite construction. ∎

### Useful exact specializations

When $a=0$, the determinant condition is vacuous: realizability is exactly the convex feasibility problem (6.3). In particular, index-two constraints alone introduce no residual nonsymmetric determinant test.

When $k=0$, the infinite part is index one, $K$ is absent, and the determinant is that of the full algebraic multiplier $Z_{00}$.

When $h=0$, all algebraic blocks disappear and the criterion is the standard-system KYP inequality with $X\succeq\epsilon F^TF$. This permits stationary zero-energy directions rather than imposing $X\succ0$.

When $f=0$, there is no finite-dynamics inequality. None of these cases needs a separate limiting argument.

For the completely controllable, completely observable class of [R3], the cited positive-real theorem provides the clean specialization: with $D+D^T\succeq0$, positive realness is sufficient as well as necessary. The present theorem does not assume those rank conditions and tests the supplied nonminimal realization itself.

## 7. Singular symmetric feedthrough and hidden components

Let $D_s=D+D^T$. If $D_s\not\succeq0$, the strict-equivalence problem is infeasible because $D_s$ is a principal block of every KYP matrix and transforms only by port congruence.

Suppose $D_s\succeq0$. Choose an orthogonal port basis $V=(V_+,V_0)$, where $V_0$ spans $\ker D_s$, and

$$
 V_+^TD_sV_+=\Sigma\succ0.
$$

Put

$$
 H_Q=-A^TQ-Q^TA,\qquad L_Q=C^T-Q^TB.
$$

Then the full KYP inequality is equivalent to

$$
 L_QV_0=0,
 \qquad
 \begin{bmatrix}
 H_Q&L_QV_+\\V_+^TL_Q^T&\Sigma
 \end{bmatrix}\succeq0.
 \qquad\text{(7.1)}
$$

**Proof.** After orthogonal congruence, the kernel-port diagonal block is zero. Lemma 4.1 makes its off-diagonal block zero as well, yielding the equality in (7.1). Removing these zero rows and columns leaves the displayed PSD block. Reattaching them proves sufficiency. ∎

If desired, the remaining PSD condition is equivalent to its Schur complement

$$
 H_Q-L_QV_+\Sigma^{-1}V_+^TL_Q^T\succeq0,
$$

but the block form (7.1) keeps the representation affine. When $D_s=0$, the condition is simply $L_Q=0$ and $H_Q\succeq0$.

Thus no inverse of a singular feedthrough is needed, and no controllability/observability recursion is assumed. Further known zero faces can be compressed in the same way. A numerical facial-reduction scheme may use these identities, but termination of the exact theorem does not depend on a successful floating-point facial-reduction implementation.

All finite and infinite states remain present in (6.2). In particular, unreachable unstable finite states, or invisible high-index constraints, cannot be dropped to manufacture a positive answer. Such components are tested through $F$, the full KYP matrix, and the infinite-index condition. Conversely, the theorem allows a singular finite metric when its kernel consists of directions annihilated by $F$, with the additional restrictions imposed by the full KYP coupling.

## 8. A finite convex-feasibility and interpolation algorithm

A determinant inequality alone is not the final algorithm. This section decides it with finitely many convex subproblems.

### 8.1 Exact encoding of strict affine inequalities

For a real affine scalar expression $\ell$,

$$
 \ell>0
 \iff\exists w\in\mathbb R:\quad
 \begin{bmatrix}\ell&1\\1&w\end{bmatrix}\succeq0.
 \qquad\text{(8.1)}
$$

Indeed, PSD gives $\ell\ge0$ and $\ell w\ge1$, hence $\ell>0$. Conversely choose $w=1/\ell$. The condition $\epsilon>0$ and all strict linear inequalities below can therefore be expressed using additional ordinary affine LMIs. No Slater condition is assumed.

### 8.2 Construct an affine-spanning family of feasible projections

Use an exact convex-feasibility oracle that returns a feasible point or a correct infeasibility decision for affine matrix inequalities. Section 9 explains a terminating real-algebraic implementation.

First decide whether $\mathscr F$ is empty. If not, obtain a feasible tuple $p_0$, with algebraic block $Z_0$. Maintain feasible tuples $p_0,\ldots,p_r$, their blocks $Z_0,\ldots,Z_r$, and

$$
 \mathcal V_r=\mathrm{span}\{\mathrm{vec}(Z_i-Z_0):1\le i\le r\}.
$$

Choose a basis $h_1,\ldots,h_b$ of $\mathcal V_r^\perp$. For each basis vector, test feasibility of the original set $\mathscr F$ with either one of the additional conditions

$$
 h_j^T\mathrm{vec}(Z_{00}-Z_0)>0,
 \quad\text{or}\quad
 h_j^T\mathrm{vec}(Z_{00}-Z_0)<0.
 \qquad\text{(8.2)}
$$

A successful query supplies a new feasible tuple whose projected difference is outside $\mathcal V_r$, so the dimension strictly increases. If every query fails, every feasible projected difference is orthogonal to every $h_j$. Therefore

$$
 \mathrm{aff}(\mathscr Z)=Z_0+\mathcal V_r.
 \qquad\text{(8.3)}
$$

The dimension can increase at most $a^2$ times. This procedure consequently terminates after finitely many exact convex-feasibility calls. A loose bound is $1+2a^2(a^2+1)$. Discard any redundant points, obtaining $p_0,\ldots,p_d$ whose projected differences form a basis, where $d\le a^2$.

Neither closedness of $\mathscr Z$ nor attainment of a linear optimum is required. We ask feasibility of a strict inequality, not whether an optimum is attained.

### 8.3 Deterministic finite grid inside a feasible simplex

For $a=0$, accept immediately after feasibility. For $a\ge1$, define

$$
 p(t_1,\ldots,t_d)
 =\det\left(Z_0+\sum_{i=1}^d t_i(Z_i-Z_0)\right).
 \qquad\text{(8.4)}
$$

This polynomial has total degree at most $a$, hence degree at most $a$ in each coordinate. When $d=0$, just evaluate $\det Z_0$.

For $d>0$, set

$$
 L=d(a+1)+1,
 \qquad t_i\in\left\{\frac1L,\frac2L,\ldots,\frac{a+1}{L}\right\}.
 \qquad\text{(8.5)}
$$

Each of the $(a+1)^d$ grid points has $t_i>0$ and $\sum_i t_i<1$. Thus

$$
 p_*=(1-\textstyle\sum_i t_i)p_0+\sum_i t_i p_i
 \qquad\text{(8.6)}
$$

is a convex combination of **complete feasible tuples**, not just of their algebraic blocks. It remains in $\mathscr F$, with a strictly positive combined $\epsilon$.

If (8.4) is nonzero at any grid point, (8.6) is a realizability certificate.

If it vanishes at every grid point, it is the zero polynomial. To see this, a univariate polynomial of degree at most $a$ vanishing at $a+1$ distinct points is zero. Apply this statement successively in each coordinate, or induct on $d$ by viewing the polynomial as a polynomial in its final variable. Its coefficient polynomials vanish on the preceding product grid and therefore vanish identically. Hence $p\equiv0$ on $\mathbb R^d$.

By (8.3), every element of $\mathscr Z$ is covered by that affine parametrization. Therefore every feasible $Z_{00}$ is singular. This proves a negative answer, not merely failure to find a nonsingular sample.

This proves Theorem 6.1(ii)⇔(iii) and gives an explicit finite determinant construction.

### 8.4 Complete algorithm

1. Verify regularity and compute finite/infinite coordinates using Section 10 or a certified deflating-subspace method.
2. If the infinite index exceeds two, return **NO: index obstruction**.
3. Build (6.2)–(6.3), optionally compressing singular-feedthrough zero rows using (7.1).
4. If $\mathscr F$ is empty, return **NO: reduced convex infeasibility**.
5. If $a=0$, use any feasible tuple. Otherwise construct an affine-spanning family by Section 8.2 and apply the grid in Section 8.3. If the polynomial vanishes identically, return **NO: every feasible index-one multiplier is singular**.
6. For a successful tuple, recover the original-coordinate metric
   $$
   Q=L^TQ_cT^{-1},
   \qquad\text{(8.7)}
   $$
   where $E_c=LET,A_c=LAT$, with no port change needed for the split.
7. Construct the pH matrices by (4.5), or use the stronger dissipative completion in Section 11.
8. Verify the exact factorization identities, energy PSD condition, and weighted dissipation PSD condition in the original coordinates.

Correctness follows from Theorem 6.1 and the proved affine-hull/interpolation argument. Termination follows from finite-dimensional coordinate operations, at most $a^2$ hull extensions, a finite grid, and the exact feasibility procedures of Section 9. Invariance follows from Section 3; although the numerical value of $\epsilon$ depends on coordinates, its existence is equivalent to an invariant kernel condition.

## 9. Exact computability, certificates, and limits

For rational or real-algebraic input, all finite coordinate operations in Section 10 are effective over the real algebraic numbers. A symmetric real matrix is PSD exactly when all its principal minors are nonnegative. Thus each finite collection of affine LMIs and affine equalities, including the lifted strict constraints (8.1), is a finite system of polynomial equalities and inequalities.

The standard decision theorem for real closed fields supplies a terminating feasibility and sample-point procedure for such systems with real-algebraic coefficients. This is the one general real-algebraic algorithmic fact used here, not a newly proved complexity result. A complete exact nonlinear-real-arithmetic procedure provides an alternative implementation; [R7] describes the relevant family of methods.

For clarity, the PSD principal-minor criterion also has a short proof. Necessity follows by restriction to principal submatrices. Conversely, if every principal minor of a symmetric $M$ is nonnegative, then every leading principal minor of $M+tI$ is strictly positive for $t>0$: expanding the determinant in $t$ gives nonnegative sums of principal minors and a positive leading term. The positive-definite Sylvester criterion yields $M+tI\succ0$; letting $t\downarrow0$ gives $M\succeq0$.

Positive outputs have finite certificates: exact entries of $Q$ and pH factors, factorization identities, and exact PSD checks. Negative outputs of the theoretical algorithm consist of an index obstruction, exact infeasibility decisions, or affine-hull containment decisions plus a determinant polynomial identity. The bundled SMT solver's `unsat` answers are reproducible computational evidence; they are not packaged as independently machine-checked unsatisfiability proof objects. The universal theorem is proved analytically above, and explicit negative examples below have direct mathematical proofs.

For arbitrary real input without a finite effective representation, the iff theorem remains valid over $\mathbb R$; executable termination is stated relative to exact arithmetic and feasibility oracles. One cannot promise a finite Turing encoding for an arbitrary black-box real number.

### Numerical interpretation

A numerical front end should use reordered generalized Schur and rank-revealing staircase methods to identify finite and infinite deflating subspaces. Rank decisions determine $h$, the infinite chain lengths, $k$, and $a$. At index at most two, the concrete tests are the reduced affine LMIs, $\lambda_{\min}(K)>0$, and $\sigma_{\min}(Z_{00})>0$, together with the finite dominance bound.

Exact nullspaces or certified enclosures are needed to turn boundary-sensitive numerical rank decisions into proofs. The strict margin $\epsilon$, singular values, transformation condition numbers, KYP residuals, and distance from nonzero eigenvalues to zero are useful diagnostics, but a small residual alone is not a certificate. Weakly infeasible SDPs and nearly singular algebraic multipliers require additional care. The power/Fitting construction in Section 10 is an exact-arithmetic reference method, not a recommended floating-point algorithm.

This document does not claim a globally well-conditioned, polynomial-time numerical method, nor an even-pencil implementation equivalent in numerical sophistication to all the algorithms in [R4]. Those are stronger algorithm-engineering and complexity questions, not missing implications in the exact iff theorem.

## 10. A real exact finite/infinite construction

This section proves that the coordinates required above can be found constructively without knowing eigenvalues in advance.

The scalar polynomial $\det(\alpha E-A)$ has degree at most $n$ and is not identically zero. Therefore one of the integers $\alpha\in\{0,1,\ldots,n\}$ makes $\alpha E-A$ invertible. Define

$$
 R_0=(\alpha E-A)^{-1}E.
 \qquad\text{(10.1)}
$$

The Fitting decomposition is

$$
 \mathbb R^n=\mathrm{im}R_0^n\oplus\ker R_0^n.
 \qquad\text{(10.2)}
$$

To verify it, ranks and nullities of powers stabilize by exponent $n$. If $w=R_0^nv$ and $R_0^nw=0$, then $R_0^{2n}v=0$; stabilization gives $v\in\ker R_0^n$, so $w=0$. Rank-nullity gives the direct sum. Both summands are invariant, the restriction to the image is invertible, and the restriction to the kernel is nilpotent.

Choose a basis matrix $S_0$ for the two spaces, so

$$
 S_0^{-1}R_0S_0=\mathrm{diag}(R_f,R_i),
$$

with $R_f$ invertible and $R_i$ nilpotent. State transformation $S_0$ and left transformation $S_0^{-1}(\alpha E-A)^{-1}$ produce

$$
 E'=\mathrm{diag}(R_f,R_i),\qquad
 A'=\mathrm{diag}(\alpha R_f-I,\alpha R_i-I).
$$

Since $\alpha R_i-I$ is invertible, a further left multiplication by

$$
 \mathrm{diag}\big(R_f^{-1},(\alpha R_i-I)^{-1}\big)
$$

gives (5.1), with

$$
 F=\alpha I-R_f^{-1},\qquad
 N=(\alpha R_i-I)^{-1}R_i.
 \qquad\text{(10.3)}
$$

The last matrix is nilpotent: its first factor is a polynomial in $R_i$, commutes with $R_i$, and is invertible. Test successive powers to determine the infinite index, or simply test $N^2=0$ before constructing the reduced problem.

If $N^2=0$, choose a basis $a_1,\ldots,a_k$ of $\mathrm{im}N$. Solve $Nb_i=a_i$. Complete the $a_i$ to a basis $(a_i,c_j)$ of $\ker N$. Then

$$
 (a_1,\ldots,a_k,b_1,\ldots,b_k,c_1,\ldots,c_a)
$$

is a basis: applying $N$ to a zero linear combination first annihilates all $b_i$ coefficients, and independence in $\ker N$ annihilates the rest. This produces (6.1). Use the inverse basis on the equation side to preserve $A_\infty=I$.

All operations are finite rank, nullspace, inverse, or linear-system operations. For rational input, all these coordinates can be chosen rational. No numerical Jordan decomposition is required.

## 11. Optional stronger unweighted dissipation

The target only requires the weighted PSD matrix in (1.5). In fact, after a suitable completion on the unused part of the range, one can arrange

$$
 \begin{bmatrix}R&P\\P^T&S\end{bmatrix}\succeq0
 \qquad\text{(11.1)}
$$

as well. This completion is consistent with the already-known singular-metric observation in [R2, Remark 12].

Define

$$
 \mathbb Q=\mathrm{diag}(Q,I_m),\quad
 \mathbb A=\begin{bmatrix}A&B\\-C&-D\end{bmatrix},
$$

$$
 \Theta_0=\mathbb A\mathbb Q^\dagger,\qquad
 \Pi=\mathbb Q\mathbb Q^\dagger,
$$

$$
 \Theta=\Theta_0-\Pi\Theta_0^T(I-\Pi).
 \qquad\text{(11.2)}
$$

Kernel compatibility gives $\Theta_0\mathbb Q=\mathbb A$. Since $(I-\Pi)\mathbb Q=0$, also $\Theta\mathbb Q=\mathbb A$.

In the orthogonal decomposition $\mathrm{im}\mathbb Q\oplus\ker\mathbb Q^T$, the identity $\Theta_0=\Theta_0\Pi$ gives

$$
 \Theta_0=\begin{bmatrix}T_{11}&0\\T_{21}&0\end{bmatrix},\qquad
 \Theta=\begin{bmatrix}T_{11}&-T_{21}^T\\T_{21}&0\end{bmatrix}.
$$

The KYP inequality is

$$
 -\mathbb Q^T(\Theta_0+\Theta_0^T)\mathbb Q\succeq0.
$$

Because $\mathbb Q$ maps onto its image, this forces $T_{11}+T_{11}^T\preceq0$. Thus $\Theta+\Theta^T\preceq0$. Write

$$
 \tfrac12(\Theta-\Theta^T)
 =\begin{bmatrix}J&G\\-G^T&-N_p\end{bmatrix},
$$

$$
 -\tfrac12(\Theta+\Theta^T)
 =\begin{bmatrix}R&P\\P^T&S\end{bmatrix}\succeq0.
$$

The identity $\Theta\mathbb Q=\mathbb A$ gives all the pH equations. This is the completion implemented by default in the reference code.

## 12. A larger behavior-preserving equivalence: output gauges

Strict equivalence does not include every coefficient change that leaves the full trajectory output unchanged. A useful additional operation is

$$
 \widetilde C=C-W^TA,\qquad
 \widetilde D=D-W^TB,\qquad E^TW=0,
 \qquad\text{(12.1)}
$$

where $W\in\mathbb R^{n\times m}$. For every trajectory,

$$
 \widetilde Cx+\widetilde Du
 =y-W^T(Ax+Bu)
 =y-W^TE\dot x=y.
 \qquad\text{(12.2)}
$$

Therefore (12.1) preserves all states, every admissible trajectory, the actual output, and the supply rate. It has inverse gauge $-W$. This operation is already used in the descriptor realization literature [R3]; it is stated explicitly rather than silently included in “equivalence.”

### Theorem 12.1: exact characterization with gauges

For the equivalence class generated by (3.1) and (12.1), Theorem 6.1 and the finite algorithm remain valid after adding the affine variable $W$, the equality $E_c^TW=0$, and replacing $(C_c,D)$ in the KYP matrix by

$$
 (C_c-W^TA_c,\ D-W^TB_c).
 \qquad\text{(12.3)}
$$

**Proof.** The state principal blocks and the energy condition used in Theorem 5.2 are unchanged. The output-kernel implications refer to the corrected output, which is exactly the one needed for factorization. All entries of the corrected KYP matrix are affine in $(Q,W)$, so the feasible set remains convex. The affine-hull and interpolation proof applies without alteration to its $Z_{00}$ projection.

Under (3.1), a gauge transforms as $\widehat W=L^{-T}WV$. Direct substitution gives $\widehat E^T\widehat W=0$ and the transformed corrected coefficients. Gauges compose by addition when the state equations are fixed; transporting them across strict equivalences shows that any finite sequence can be represented by a strict equivalence and one gauge. Conversely a feasible $W$ and metric explicitly construct an allowed gauged pH realization. ∎

This is not a claim to characterize every imaginable equivalence of descriptor behaviors. State/input mixing, dynamic port transformations, and changes in state dimension remain outside the declared classes.

## 13. Explicit positive and negative certificates

### 13.1 Smallest forced-singular algebraic obstruction

Take the scalar system

$$
 E=0,\quad A=B=1,\quad C=D=0.
$$

The pencil is regular and has index one. For a scalar $q$,

$$
 \mathcal W(q)=\begin{bmatrix}-2q&-q\\-q&0\end{bmatrix}.
$$

PSD forces its off-diagonal entry to vanish, so $q=0$. A KYP solution exists, but $A\ker Q\ne0$. Equivalently, the only feasible algebraic multiplier is singular. Hence there is no strict-equivalence pH realization.

Nevertheless the trajectory equations are $x=-u,y=0$; the system is passive with zero storage and has zero transfer function. This demonstrates why neither passivity nor a positive-real transfer function suffices for the fixed nonminimal realization. Allowing the output gauge changes the answer: taking $W=-1$ gives $\widetilde C=\widetilde D=1$, and $q=-1$ is a valid pH metric. Thus the declared equivalence class materially affects the answer.

### 13.2 A nonsingular skew algebraic multiplier must not be rejected

Let

$$
 J_0=\begin{bmatrix}0&1\\-1&0\end{bmatrix},\qquad
 E=0_{2\times2},\quad A=B=I_2,\quad C=J_0,\quad D=0.
$$

Since $D+D^T=0$, the KYP off-diagonal block must vanish, forcing $Q=C=J_0$. This matrix is invertible and $\mathcal W(Q)=0$, so pH realization exists. However $-Q-Q^T=0$. Replacing the determinant test by strict negativity of the symmetric part would give a false negative.

In odd dimension, the same construction with a skew $C$ forces a singular skew $Q=C$, and therefore fails the compatibility condition. This also illustrates why a common-kernel shortcut for arbitrary nonsymmetric affine families is unsafe.

### 13.3 Lossless index-two differentiator

Let

$$
 H_0=\begin{bmatrix}0&1\\0&0\end{bmatrix},\qquad
 E=H_0,\quad A=I_2,\quad B=\begin{bmatrix}0\\1\end{bmatrix},
$$

$$
 C=\begin{bmatrix}-1&0\end{bmatrix},\quad D=0,\quad Q=J_0.
$$

Then

$$
 E^TQ=\mathrm{diag}(0,1),\qquad \mathcal W(Q)=0.
$$

A pH realization is $J=-J_0,R=0,G=B,P=0,S=N_p=0$. The equations give $y=\dot u$ for smooth inputs. This is an admissible index-two realization, not a proper standard-state-space system.

### 13.4 One explicit example combines nonminimality, index two, singular feedthrough, and singular Q

Use $H_0,J_0$ above and $X_0=\mathrm{diag}(0,1)$. Define

$$
 E=\mathrm{diag}(I_2,H_0,0),\quad
 A=\mathrm{diag}(H_0,I_2,1),\quad
 B=e_4,\quad C=-e_3^T,\quad D=0.
 \qquad\text{(13.1)}
$$

An exact certificate is

$$
 Q=\mathrm{diag}(X_0,J_0,-1),
$$

$$
 J=\mathrm{diag}(J_0,-J_0,0),\quad
 R=\mathrm{diag}(0,0,0,0,1),
$$

$$
 G=e_4,\quad P=0,\quad S=N_p=0.
 \qquad\text{(13.2)}
$$

Direct multiplication gives

$$
 (J-R)Q=A,\qquad (G+P)^TQ=C,
$$

$$
 E^TQ=\mathrm{diag}(0,1,0,1,0)\succeq0,
$$

$$
 \begin{bmatrix}Q^TRQ&Q^TP\\P^TQ&S\end{bmatrix}
 =\mathrm{diag}(0,0,0,0,1,0)\succeq0.
 \qquad\text{(13.3)}
$$

Here $Q$ has rank four rather than five; $D+D^T=0$; the finite two-state nilpotent subsystem is unactuated and unobserved; and both index-two and index-one infinite blocks are present. No state was removed to obtain the certificate.

### 13.5 Determinant interpolation can need a convex combination

Take feasible projected matrices $Z_0=\mathrm{diag}(1,0)$ and $Z_1=\mathrm{diag}(0,1)$. Both vertices are singular, but their strictly interior combinations are nonsingular. With $a=2,d=1$, the prescribed grid first tests $t=1/4$, giving determinant $3/16$. Testing only the affine-spanning vertices would not be a correct algorithm.

## 14. Reproducibility and exact test results

The public repository contains:

- `code/phdae.py`: exact finite/infinite decomposition, reduced symbolic KYP problem, pH factor construction, exact PSD and identity checks, direct polynomial feasibility encoding, and determinant-grid evaluation for a supplied affine-spanning family.
- `code/z3_capi.py`: a minimal interface to the public libz3 C API; no Python `z3-solver` package is required when the shared library is available.
- `code/run_tests.py`: exact certificate and adversarial regression suite.
- `code/verify_certificates.py`: a separate rational certificate checker requiring only SymPy; it imports neither the solver interface nor `phdae.py`.
- `checks/`: saved SMT queries, solver outputs, system data, and positive exact certificates.
- `requirements.txt`: the pinned Python dependency for the verification scripts.

From the repository root, the independent certificate checker can be run with

```bash
python -m pip install -r requirements.txt
python code/verify_certificates.py
```

The larger suite is run with

```bash
python code/run_tests.py
```

The reference implementation uses SymPy and, for the larger suite, a system installation of the Z3 shared library (`libz3`). `code/z3_capi.py` searches for the library automatically; `Z3_LIBRARY` can be set to its full path if necessary.

**Implementation boundary:** `solve_reduced` implements the reduced criterion directly as exact polynomial feasibility using all PSD principal minors and $\det Z_{00}\ne0$. This is an alternative exact decision route to Section 8, not an implementation of its affine-hull oracle loop. `simplex_determinant_search` implements the interpolation grid but assumes its supplied feasible samples have the stated affine-spanning property; it does not certify that assumption. The full convex-oracle procedure is specified and proved in Section 8. The general gauge theorem is proved in Section 12; the reference solver does not automatically search all gauges.

With resource limits enabled, a solver timeout or unresolved query is returned as **UNKNOWN**, never as infeasibility. Setting `timeout_ms=None` removes the explicit solver timeout, but finite physical resources still limit a real program. Positive returned models are substituted into the original matrices and independently checked with exact arithmetic. An unresolved algebraic sign raises an error instead of being accepted on numerical evidence.

The observed suite completed **20 test groups**, including eleven named feasible/infeasible systems, the explicit mixed certificate, strict-equivalence tests, twelve rationally scrambled decompositions, three full-Q cross-block counterexample searches, eight reduced-versus-unreduced feasibility comparisons, determinant interpolation checks, and an output-gauge scope example. All groups passed. A separate checker also verified all five saved rational positive certificates against their saved original system matrices, without invoking the solver or importing `phdae.py`. These finite tests support implementation checking; the universal iff theorem does not follow from test counts.

## 15. Verification checklist

- [x] The original pencil is assumed regular, not necessarily minimal or index one.
- [x] The target pH convention is the one with singular Q allowed and weighted dissipation PSD.
- [x] Admitted state/equation/port transformations and their inverses are explicit.
- [x] Full trajectory equivalence and power preservation are proved.
- [x] No controllability/observability reduction deletes hidden states.
- [x] The known KYP-plus-kernel theorem is credited and independently proved.
- [x] The nilpotent Lyapunov lemma has a complete proof.
- [x] Vanishing of the finite/infinite upper block holds for every KYP solution.
- [x] The off-diagonal block annihilates ker X; this is used in the kernel proof.
- [x] Infinite-block invertibility is necessary, not silently assumed.
- [x] Index greater than two is rejected by an exact proved obstruction.
- [x] Singular finite energy is handled by X ≥ epsilon F^T F, not by an invalid Cholesky inverse.
- [x] The index-two block form and determinant identity are proved, including singular intermediate K.
- [x] Singular symmetric feedthrough is treated without an inverse of D+D^T.
- [x] The remaining determinant condition has a finite decision procedure.
- [x] The affine-hull procedure handles nonclosed/unbounded projections and strict inequalities.
- [x] The finite interpolation grid lies inside a feasible simplex of full tuples.
- [x] Exact algorithm termination is stated with an appropriate input/arithmetic model.
- [x] Explicit pH matrices are constructed and verified in original coordinates.
- [x] The optional stronger unweighted dissipative completion is proved.
- [x] Output gauges are separated from strict equivalence and have their own iff theorem.
- [x] Exact scripts and observed outputs are preserved; UNKNOWN is not treated as NO.
- [x] Theoretical oracle algorithm and implemented polynomial solver are not conflated.
- [x] Earlier descriptor realization results and the 2025 publication update are acknowledged.
- [x] No historical-priority, polynomial-complexity, uniform-stability, or independent-review claim is made.

## 16. References

**[R1]** C. Beattie, V. Mehrmann, H. Xu, H. Zwart. *Linear port-Hamiltonian descriptor systems*. Mathematics of Control, Signals, and Systems 30, Article 17 (2018). DOI: **10.1007/s00498-018-0223-3**. Preprint: arXiv:1705.09081.

**[R2]** K. Cherifi, H. Gernandt, D. Hinsen. *The difference between port-Hamiltonian, passive and positive real descriptor systems*. Mathematics of Control, Signals, and Systems 36 (2024), 451–482; first published online in 2023. DOI: **10.1007/s00498-023-00373-2**. Particularly Proposition 11, Remark 12, and Section 4.

**[R3]** D. Chu, V. Mehrmann. *Port-Hamiltonian representations of positive real descriptor systems*. Automatica 180 (2025), 112456. DOI: **10.1016/j.automatica.2025.112456**. Preprint under the title *Port-Hamiltonian realizations of positive real descriptor systems*, arXiv:2408.14115. Particularly Theorem 7, Corollary 8, Theorem 10, and Corollary 11.

**[R4]** C. Beattie, V. Mehrmann, H. Xu. *Port-Hamiltonian realizations of non-minimal linear time-invariant systems*. Mathematics of Control, Signals, and Systems 38 (2026), 1–41. DOI: **10.1007/s00498-025-00432-w**. Preprint: arXiv:2201.05355. Particularly Remarks 30–31.

**[R5]** C. Mehl, V. Mehrmann, M. Wojtylak. *Linear algebra properties of dissipative Hamiltonian descriptor systems*. SIAM Journal on Matrix Analysis and Applications 39 (2018), 1489–1519. DOI: **10.1137/18M1164275**. Preprint: arXiv:1801.02214.

**[R6]** D. Chu, V. Mehrmann. *Regularization of port-Hamiltonian descriptor systems*. Electronic Transactions on Numerical Analysis 65 (2026), 254–270. Published online June 26, 2026. DOI: **10.1553/etna_vol65s254**. Preprint: arXiv:2509.02715.

**[R7]** D. Jovanović, L. de Moura. *Solving Non-Linear Arithmetic*. Automated Reasoning, IJCAR 2012, Lecture Notes in Computer Science 7364. The exact reference implementation uses the public Z3 C API and the `qfnra-nlsat` tactic. The finite mathematical decision claim is the standard real-closed-field polynomial feasibility theorem, not a performance guarantee for a bounded solver run.
