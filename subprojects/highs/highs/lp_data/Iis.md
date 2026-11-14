# HiGHS irreducible infeasibility system (IIS) facility

Further to documentation in https://ergo-code.github.io/HiGHS/stable/guide/advanced/

Source of IIS techniques is

J. W. Chinneck, "Feasibility and Infeasibility in Optimization". Algorithms and Computational Methods, Springer, New York, 2008. [DOI 10.1007/978-0-387-74932-7](https://doi.org/10.1007/978-0-387-74932-7)

The IIS search is rooted in `Highs::getIisInterface()`, which first
checks whether the `Highs::model_status_` is
`HighsModelStatus::kOptimal` or `HighsModelStatus::kUnbounded`, in
which case the model is feasible so no IIS exists. Otherwise, the
trivial check (for inconsistent bounds or empty infeasible rows) -
performed by `HighsIis::trivial()` - and infeasible rows based on row
value bounds - performed by `HighsIis::rowValueBounds()` - are
performed. If `Highs::options_.iis_strategy` is `kIisStrategyLight`
then `Highs::getIisInterface()` returns.

The "full" IIS calculation operates in two phases: after a set of
mutually infeasible rows has been identified, this is reduced to an
IIS. 

## Finding a set of mutually infeasible rows

The set of mutually infeasible rows can be found in two ways. For
a well-documented example of an IIS and its calculation, see unit test
`lp-get-iis-galenet` in `check/TestIis.cpp`.

### Using a dual ray

If it is known that the model is infeasible, then the simplex
solver may have identified a dual ray. If there is a dual ray then its
nonzeros correspond to a set of mutually infeasible constraints. If
there is no dual ray - as might happen if the model's infeasibility
has been identified in presolve - then the incumbent model is solved
with `Highs::options_.presolve` = "off". Unfortunately the "ray route"
is not robust, so currently switched off.

### Using the elasticity filter

Currently, the only route to finding a set of mutually infeasible rows
is to perform an elasticity filter calculation [Chinneck
p101-104]. This is done in `HighsIis::elasticityFilter`. This method
is more general than is necessary for finding the set of mutually
infeasible rows for an IIS calculation, and can be called directly
using `Highs::feasibilityRelaxation`.

The essence of the `HighsIis::elasticityFilter` is that it allows
lower bounds, upper bounds and RHS values to be violated. There are
penalties for doing so that can be global for each of these three
cases, or local to each column lower bound, upper bound and row
bound. The "elasticity LP" is constructed by adding non-negative
elastic variables to transform the constraints from
$$
L <= Ax <= U;\qquad l <= x <= u
$$
to
$$
L <= Ax + e_L - e_U <= U;\qquad l <=  x + e_l - e_u <= u,
$$
where the elastic variables are not used if the corresponding bound is
infinite or the local/global penalty is negative. If a bound on $x$
can be violated, then it is removed from $x$. This is because the
violation is modelled using a row $l <= x + e_l - e_u <= u$ containing
at least one of $e_l$ and $e_u$, depending on whether both bounds can
be violated. The objective of the elasticity LP is the linear function
of the elastic variables given by the local/global penalties. Note
that the model modifications required to achieve this formulation, are
made by calls to methods in the `Highs` class so that the value of any
initial basis is maintained.

For the purposes of IIS calculation, since an infeasible subset of
rows is required, only row bounds are allowed to be violated, and the
penalty in each case is unity. Hence `HighsIis::elasticityFilter` uses
`global_lower_penalty=-1`, `global_upper_penalty=-1` and
`global_rhs_penalty=1`. Each finite row bound has a corresponding
elastic variable.

When the elasticity LP is solved, if all elastic variables have an
optimal value of zero, then the original LP is feasible. Clearly the
algorithm to find an infeasible subset of rows then
terminates. Otherwise, the aim of the algorithm is to identify a set
of elastic variables that, if fixed to zero, makes the elasticity LP
infeasible. The algorithm proceeds as follows. After solving the
elasticity LP, each elastic variable whose optimal value is positive
is fixed at zero, and the elasticity LP is re-solved. This process
continues until the elasticity LP is infeasible. Thus, the elastic
variables that have been fixed at zero correspond to an infeasible
subset of rows. Note that the algorithm must continue until the
elasticity LP is infeasible as, for example, the initial set of
elastic variables that have positive values may not correspond to an
infeasible subset of rows. It is possible that these constraints can
be satisfied, but only at the cost of a greater total violation of
some other constraint(s) that were satisfied.

## Reducing to an IIS

Given an infeasible subset of rows, the method `HighsIis::getData`
forms the LP corresponding to the subset of rows and the columns with
nonzeros in those rows, and the row bounds corresponding to the fixed
elastic variables. This LP is then passed to `HighsIis::compute`,
which eliminates rows and columns from the LP until this is no longer
possible, so the subset of rows and the columns is
irreducible. Depending on whether the user desires a minimum number of
rows or columns in the final IIS, it is possible to prioritise the
removal of rows or columns according to the value of
`Highs::options_.iis_strategy`. With the aim of finding a minimum
number of rows in the final IIS, removal technique considers each
finite row bound in turn. If it is relaxed and the the LP is still
infeasible, then the bound can be removed. If this leads to both
bounds on a row being infinite, then the row itself can be removed. A
similar pass through the columns is performed before or after the row
pass, depending on whether removal of columns or rows has priority.
