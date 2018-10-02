/**
 * @name Pythagorean calculation with sub-optimal numerics
 * @description Calculating the length of the hypotenuse using the standard formula may lead to overflow.
 * @kind problem
 * @tags accuracy
 * @problem.severity warning
 * @sub-severity low
 * @precision medium
 * @id py/pythagorean
 */

import python

predicate squareOp(BinaryExpr e) {
  e.getOp() instanceof Pow and e.getRight().(IntegerLiteral).getN() = "2"
}

predicate squareMul(BinaryExpr e) {
  e.getOp() instanceof Mult and e.getRight().(Name).getId() = e.getLeft().(Name).getId()
}

predicate squareRef(Name e) {
  e.isUse() and
  exists(SsaVariable v, Expr s |
    v.getVariable() = e.getVariable() |
    s = v.getDefinition().getNode().getParentNode().(AssignStmt).getValue() and
    square(s)
  )
}

predicate square(Expr e) {
  squareOp(e)
  or
  squareMul(e)
  or
  squareRef(e)
}

from
  Call c,
  BinaryExpr s
where
  c.getFunc().toString() = "sqrt" and
  c.getArg(0) = s and
  s.getOp() instanceof Add and
  square(s.getLeft()) and square(s.getRight())
select
  c, "Pythagorean calculation with sub-optimal numerics"