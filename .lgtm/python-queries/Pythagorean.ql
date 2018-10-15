/**
 * @name Pythagorean calculation with sub-optimal numerics
 * @description Calculating the length of the hypotenuse using the standard formula may lead to overflow.
 * @kind problem
 * @tags accuracy
 * @problem.severity warning
 * @sub-severity low
 * @precision high
 * @id py/pythagorean
 */

import python

// Three different ways to write squares:
// a**2
predicate squareOp(BinaryExpr e) {
  e.getOp() instanceof Pow and e.getRight().(IntegerLiteral).getN() = "2"
}

// a*a
predicate squareMul(BinaryExpr e) {
  e.getOp() instanceof Mult and e.getRight().(Name).getId() = e.getLeft().(Name).getId()
}

// v, where v is assigned to a square
predicate squareRef(Name e) {
  e.isUse() and
  exists(SsaVariable v, Expr s |
    v.getVariable() = e.getVariable() |
    s = v.getDefinition().getNode().getParentNode().(AssignStmt).getValue() and
    square(s)
  )
}

// a square is either of the three
predicate square(Expr e) {
  squareOp(e)
  or
  squareMul(e)
  or
  squareRef(e)
}

// a hypot is sqrt on a sum of two squares
predicate hypot(Call c) {
  c.getFunc().toString() = "sqrt" and
  exists(BinaryExpr s |
    c.getArg(0) = s and
    s.getOp() instanceof Add and
    square(s.getLeft()) and square(s.getRight())
  )
}

from Call c
where
  hypot(c)
select
  c, "Pythagorean calculation with sub-optimal numerics"
