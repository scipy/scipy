# `Highs::run()`

`Highs::run()` has evolved a great deal since it was first created to
"solve" the LP in the `HighsLp` instance `Highs::lp_`. As well as
solving more problem classes, and using more solvers, features have
been added inelegantly. The simplicity of having a single
`Highs::run()` method, with different features obscured within it will
be replaced by multiple, nested, "solve" methods which are explicit in
their features and very much simpler.

The only refactoring came when the multiple objective code was added:
`Highs::optimizeModel()` inherited the content of `Highs::run()` so
that a single call to`Highs::run()` could perform multiple
optimizations.

Other developments that have been implemented inelegantly are

### Actioning executable run-time options via `Highs::run()`

As [#2269](https://github.com/ERGO-Code/HiGHS/issues/2269)
highlighted, users of `cvxpy` can only execute `Highs::run()`, so the
following actions that were previously in `app/RunHighs.cpp`, are now
performed in `Highs::run()`

- Read from a solution and/or basis file
- Write out the model
- Write out the IIS model
- Write out the solution and/or basis file. 

There is still one action in`app/RunHighs.cpp` that should be performed in `Highs::run()`

- Write out the presolved model

These "HiGHS files" actions must only be performed at the "top level"
of `Highs::run()`, and this is achieved by caching the file options in
the `Highs` class and clearing them from options_ so that they aren't
applied at lower level calls to `Highs::run()`. They are then restored
before returning from `Highs::run()`.

### Performing user scaling

User objective and/or bound scaling is performed before assessing
whether there is excessive problem data and suggesting user objective
and bound scaling. These user scaling actions must only be performed
at the "top level" of `Highs::run()`, and this is achieved by caching
the user scaling options in the `Highs` class and clearing them from
options_ so that they aren't applied at lower level calls to
`Highs::run()`. If user scaling has been applied in a call to
`Highs::run()`, it is unapplied and the option values restored before
returning from `Highs::run()`.

### Applying "mods"

The `HighsLp` class contains data values and structures that cannot be handled explicitly by the solvers.

- If a variable has an excessively large objective cost, this is
  interpreted as being infinite, and handled in
  `Highs::handleInfCost()` by fixing the variable at its lower or
  upper bound (when finite) according to the sign of the cost and the
  sense of the optimization, and zeroing the cost. After solving the
  problem, the cost and bounds must be restored.

- If a variable is of type `HighsVarType::kSemiContinuous` or
  `HighsVarType::kSemiInteger`, it is assessed in
  `assessSemiVariables` (in `HighsLpUtils.cpp`) before reformulation
  in `withoutSemiVariables` (in `HighsLpUtils.cpp`)

  - If it is not strictly "semi" it is set to`HighsVarType::kContinuous` or `HighsVarType::kInteger`
  - If its lower bound is not positive, it is deemed to be illegal
  - If its upper bound is larger than `kMaxSemiVariableUpper` then,
    depending on the lower bound, if it is possible to reformulate it the
    upper bound is set to `kMaxSemiVariableUpper` (it is said to be "tightened". 
    Otherwise, it is deemed to be illegal

- If the relaxation is being solved, the lower bound of any
  semi-variable is set to zero in `relaxSemiVariables` (in
  `HighsLpUtils.cpp`)

All these modifications are applied in
`Highs::optimizeModelApplyMods()`, and recorded in `HighsLp::mods_`.
They are undone by `Highs::optimizeModelUndoMods()`, which is called
from `Highs::returnFromOptimizeModel()` or, when returning before
then (including after an exception), directly. The mods are applied
for each call to `Highs::optimizeModel()`, rather than once in
`Highs::optimizeHighs()`, since handling infinite costs depends on the
objective, which changes between the optimizations when there are
multiple objectives. A nested call to `Highs::optimizeModel()` from
`Highs::completeSolutionFromDiscreteAssignment()` sees the model with
mods applied, so doesn't apply or undo them.

With the plan to allow indicator constraints and SOS as generalised
disjunctive forms that will be reformulated, the handling of "mods"
needs to be refactored!

## Refactoring

The inelegance of `Highs::run()` (and `Highs::optimizeModel()`) was
exposed by
[\#2635](https://github.com/ERGO-Code/HiGHS/issues/2635). Both methods
need to be refactored. Firstly, `Highs::run()` must be refactored into
the following set of nested methods. By calling the appropriate
method, there is no need to "hide" option settings by caching and then
clearing their value. 

Refactoring `Highs::optimizeModel()` is trickier. There needs to be a
method where any "mods" are made, so that at the level below the
problem defined by the `HighsModel` class (without semi-variables) is
solved.

### `Highs::run`

This "outer" layer can contain the "HiGHS files" actions that were
previously in `app/RunHighs.cpp`, and user scaling

### `Highs::optimizeHighs()`

The next layer calls `Highs::optimizeModel()` or
`Highs::multiobjectiveSolve()` if there are multiple objectives.

### `Highs::optimizeModel()`

The next layer applies any "mods" to the `HighsModel` class (see
above) in `Highs::optimizeModelApplyMods()`, before optimizing it, and
undoes them before returning

## Observations

- Refactoring `Highs::optimizeModel()` is tricky, so is is temporarily
  renamed `Highs::calledOptimizeModel()`, and `Highs::optimizeModel()`
  is a temporary intermediate method to facilitate this is.

- The most obvious place where `Highs::run()` was called at a "lower
  level" is in the MIP solver. Since there are no "upper level"
  actions to be performed, it can call `Highs::optimizeModel()`. To
  emphasise that just an LP is being solved, `Highs::optimizeLp()` has
  been created. This is currently a call to `Highs::optimizeModel()`

## ToDo

- Move the code to write out the presolved model from `app/RunHighs.cpp` to `Highs::run()`

- Make the bound changes in
  `Highs::processedBoundsFindsInfeasibility()` a "mod" that is undone,
  rather than a permanent change to the incumbent model

- For problems with multiple objectives

  - Apply user scaling to the multiple objectives 

  - Assess the multiple objectives for extreme values

  - Accumulate subsystem solve times

- In IIS calculations, accumulate subsystem solve times




