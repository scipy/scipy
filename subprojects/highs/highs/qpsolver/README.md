# HiGHS Active set QP solver (QUASS)

This document is not a full accout of how active set methods solver
QPs (see, for example, Fletcher's _Practical Methods of Optimization_,
Chapter 10), rather some background information to help understand the
essence of the HiGHS solver. For a

## Active set vs simplex for LP

To document the active set QP solver, it is convenient to start with
observations on the active set method (ASM) for LP. Although the ASM
is equivalent to the primal simplex algorithm, this is obscured by
differences in terminology, description and underlying linear algebra.

Whilst the simplex algorithm considers an LP in standard form with
$n+m$ variables ($n$ original and $m$ slack), the ASM considers an LP
with $m+n$ constraints ($m$ original and $n$ corresponding to variable
bounds). Hence, whilst the simplex algorithm definition relates to
"variables" (that are original or slack) the ASM relates to
"constraints" (that are original or variable bounds).

Whereas the simplex method is based on the factorization of a matrix
of dimension $m$ (whose columns correspond to the basic variables) the
ASM is based on the factorization of a matrix of dimension $n$. The
rows of this matrix correspond to a set of active constraints
(original constraints at a bound and variables at a bound). As such,
the set of active constraints in the ASM is identical to the set of
nonbasic variables in the simplex algorithm. Confusingly, and due to
the role of this matrix, the set of $n$ active constraints in the ASM
is referred to as the "basis". The remaining $m$ constraints in the
ASM are analogous to the set of basic variables in the simplex
algorithm, in that their values are solved for.

For the simplex method, below is an illustration of the partition of
the original and slack variables into basic and nonbasic variables

| Basic original | Basic slack | Nonbasic slack | Nonbasic original |
| --- | --- | --- | --- |
| $B_0$ | $0$ | $I_N$ | $N_0$ |
| $B_1$ | $I_B$ | $0$ | $N_1$ |

Hence the simplex method factorizes the matrix
$\left[\begin{matrix}B_0 & 0\\ B_1 &
I_B\end{matrix}\right]\in\R^{m\times m}$.  For the ASM, below is an
illustration of the partition of the original constraints and variable
bounds into active and inactive constraints

ASM | Simplex basic | Simplex nonbasic |
| --- | --- | --- |
Active original constraints | $B_0$ | $N_0$ |
Active variable bounds| $0$ | $I_N$ |
Inactive variable bounds | $I_B$ | $0$ |
Inactive original constraints | $B_1$ | $N_1$ |

Hence the ASM factorizes the matrix $\left[\begin{matrix}B_0 & N_0\\ 0
& I_N\end{matrix}\right]\in\R^{n\times n}$.  Note that the
factorization requirement for simplex and ASM is essentially the same,
since it is restricted to the "kernel" matrix $B_0$. This is due to
the matrices $I_B$ for simplex and $I_N$ for ASM yielding singleton
pivots in Gaussian elimination.

## An active set method for QP

The QP problem solved by HiGHS is

$$
\min \quad q(x) = \dfrac{1}{2}x^TQx + c^Tx \qquad \textrm{s.t.}~ \quad L \leq Ax \leq U; \quad l \leq x \leq u,\qquad(1)
$$

where $Q$ must be positive semi-definite. Unlike LP where (if an
optimal solution exists) there is an optimal solution at a vertex of
the feasible region, in general the optimal solution of a QP is not a
vertex. This is trivially so in the case of unconstrained convex QP
problems. When a QP has linear constraints, so long as the no solution
of the unconstrained QP is feasible, there will be at least one active
constraint at an optimal solution. The null space orthogonal to these
active constraints is fundamental to the properties of an optimal
solution, and the characterisation of constraint null spaces is
critical to finding an optimal solution.

### Karush-Kuhn-Tucker (KKT) conditions for convex QP

A point $x$ is an optimal solution of (1) iff $x$ is feasible and
there is a set $\cal{A}$ of indices of active constraints (and
bounds), and Lagrange multipliers $\lambda$, such that
$Qx+c=A_{\cal{A}}^T\lambda$ (where $A_{\cal{A}}$ is the matrix
corresponding to the active constraints) and the signs of the
components of $\lambda$ are consistent with the bounds at which the
constraint is active. Note that $\cal{A}$ need only include the
indices of active constraints with nonzero Lagrange multiplier, so it
can be assumed that $A_{\cal{A}}$ is of full rank.

The aim of a QP active set method is to find $\cal{A}$ such that
$A_{\cal{A}}$ is of full rank and the KKT conditions hold. 

### Some linear algebra and calculus

For a set $\cal{A}$ and matrix $A_{\cal{A}}$, there is a "null space"
normal to the row space of $A_{\cal{A}}$. This null space plays a
major role in defining an active set method for QP. Suppose that
$r=|\cal{A}|$, so $A_{\cal{A}}\in\R^{r\times n}$, and the matrix
$Y^T\in\R^{r\times n}$ spans the row space of $A_{\cal{A}}$ such that
$A_{\cal{A}}Y=I$ (so $Y$ is the right inverse of $A_{\cal{A}}$). Let
$Z^T\in\R^{(n-r)\times n}$ span the null space of $A_{\cal{A}}$, so
$A_{\cal{A}}Z=0$.

With these definitions, note that the point $x=Yb_{\cal{A}} + Zy$,
where $b_{\cal{A}}$ are the active bounds for the constraints indexed
by $\cal{A}$, satisfies
$$
A_{\cal{A}}x=A_{\cal{A}}Yb_{\cal{A}} + A_{\cal{A}}Zy = b_{\cal{A}}.
$$
Hence $y$ parameterises $x$ in the null space.

Active set methods form and operate with $Y$ and $Z$ by identifying a
matrix $V\in\R^{(n-r)\times n}$ such that
$\displaystyle B_{\cal{A}}=\left[\begin{matrix}A_{\cal{A}}\\V\end{matrix}\right]$
is well conditioned, so that
$\left[\begin{matrix}Y&Z\end{matrix}\right]=B_{\cal{A}}^{-1}$. Note
that $Y$ and $Z$ are not formed explicitly since, for example,
$Yb_{\cal{A}}$ can be formed by factorizing $B$ and solving
$$
B_{\cal{A}}y=\left[\begin{matrix}b_{\cal{A}}\\0\end{matrix}\right]\qquad(2)
$$

Note that substituting $x=Yb_{\cal{A}} + Zy$ into the objective yields
$$
q(y) = \frac{1}{2} y^TZ^TQZy +(c+QYb_{\cal{A}})^TZy +\left(c+\frac{1}{2}QYb_{\cal{A}}\right)^TYb_{\cal{A}}
$$
Hence the minimizer of $q(y)$ is a solution of
$$
Z^TQZy = -Z^T(c+QYb_{\cal{A}})\qquad(3)
$$
The matrix $Z^TQZ$ is known as the **reduced Hessian**.

### Outline ASM algorithm

The following algorithmic definition is from Fletcher's _Practical
Methods of Optimization_, Chapter 10, which also illustrates it
qualitatively. Key to the definition is the **equality problem**, 
$$
\min~ \frac{1}{2}\delta^TQ\delta + \delta^Tg^{(k)}\quad \textrm{s.t.}~ a_i^T\delta=0, i\in\cal{A},
$$
where $g^{(k)}=Qx^{(k)}+c$. This problem minimizes the QP objective in the null space of the active constraints at a point $x^{(k)}$. Also, the set of inequalities in the QP is denoted $\cal{I}$.

1. Given feasible $x^{(1)}$ and $\cal{A}$, set $k=1$
2. If $\delta=0$ does not solve the equality proble, go to 4.
3. Compute Lagrange multipliers $\lambda^{(k)}$ for the equality
problem and use the following to determine $q$. If
$\textrm{sgn}(q)\lambda^{(k)}_q\ge0$, then terminate with
$x^*=x^{(k)}$, otherwise remove $q$ from $\cal{A}$.
$$
q=\argmin_{\cal{A}\cap\cal{I}} \textrm{sgn}(i)\lambda^{(k)}_i\quad\textrm{where}~\textrm{sgn}(i) = \begin{cases}\phantom{-}1,~~a_i^Tx=L_i\\ -1,~~a_i^Tx=U_i\end{cases}
$$
4. Solve the equality problem for $s^{(k)}$
5. Find $\alpha^{(k)}$ to solve the following **linesearch problem** and set $x^{(k+1)} = x^{(k)} + \alpha^{(k)}s^{(k)}$.
$$
\alpha^{(k)} = \min\left(1, 
\min_{i: i\not\in\cal{A}, a_i^Ts^{(k)}<0} \frac{L_i-a_i^Tx^{(k)}}{a_i^Ts^{(k)}},
\min_{i: i\not\in\cal{A}, a_i^Ts^{(k)}>0} \frac{U_i-a_i^Tx^{(k)}}{a_i^Ts^{(k)}}
\right)
$$
6. If $\alpha^{(k)}<1$ then add $p$ to $\cal{A}$, where $p\not\in\cal{A}$ yields $\alpha^{(k)}$
7. Set $k=k+1$ and go to 2.

In essence, the algorithm solves the equality problem repeatedly,
allowing constraints to become inactive if they will improve the
objective (deducing optimality if none will do so) and forcing
constraints to be active if they would be violated at the minimizer of
the equality problem. It can be seen to be a generalisation of the
primal simplex method, where nonbasic variables are allowed to move
from bounds at vertices to improve the objective, to be replaced by
basic variables which reach bounds. Since the primal simplex method is
always at a vertex, the null space is trivial, and since there is no
local minimizer along an edge of the polytope, a limiting basic
variables is always found (for an LP that is not unbounded).

### Practicalities

The algorithm assumes that $x^{(1)}$ is feasible. This can be found by
using the simplex algorithm to solve the LP feasibility problem. This
yields a vertex solution, so the initial null space is empty. The
simplex nonbasic variables yield $\cal{A}$, for which the matrix
$B_{\cal{A}}$ is nonsingular since it has the same kernel as the
simplex basis matrix.
