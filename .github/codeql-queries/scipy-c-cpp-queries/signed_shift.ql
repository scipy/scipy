/**
 * @name Signed shift
 * @description Shifting a negative number is undefined behavior,
 *              so it is risky to shift a signed number.
 * @kind problem
 * @problem.severity warning
 * @id cpp/signed-shift
 * @tags security
 *       external/cwe/cwe-758
 */

// See the "Bitwise shift operators" section here:
// https://en.cppreference.com/w/cpp/language/operator_arithmetic
import cpp
import semmle.code.cpp.rangeanalysis.SimpleRangeAnalysis

from BinaryBitwiseOperation shift, Expr lhs
where
  (shift instanceof LShiftExpr or shift instanceof RShiftExpr) and
  lhs = shift.getLeftOperand().getFullyConverted() and
  lowerBound(lhs) < 0
select shift,
  "This signed shift could cause undefined behavior if the value is negative. Type of lhs: " +
    lhs.getType().toString()
