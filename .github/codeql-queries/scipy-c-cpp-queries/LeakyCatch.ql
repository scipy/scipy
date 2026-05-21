/**
 * @name Leaky catch
 * @description If an exception is allocated on the heap, then it should be deleted when caught.
 * @kind problem
 * @problem.severity warning
 * @precision high
 * @id cpp/catch-missing-free
 * @tags efficiency
 *       correctness
 *       exceptions
 *       external/cwe/cwe-401
 */

import cpp

predicate doesRethrow(Function f) {
  exists(ReThrowExpr e | e.getEnclosingFunction() = f |
    not e.getEnclosingStmt().getParent*() instanceof CatchBlock
  )
  or
  exists(FunctionCall fc | fc.getEnclosingFunction() = f | doesRethrow(fc.getTarget()))
}

predicate deletesException(Expr expr, Parameter exception) {
  expr.getEnclosingBlock().getParent*().(CatchBlock).getParameter() = exception and
  (
    exists(FunctionCall fc | fc = expr |
      // Calling a delete function on the exception will free it (MFC's CException has a Delete function).
      fc.getQualifier() = exception.getAnAccess() and
      fc.getTarget().getName().toLowerCase().matches("%delete%")
      or
      // Passing the exception to a function might free it.
      fc.getAnArgument() = exception.getAnAccess()
      or
      // Calling a function which rethrows the current exception might cause the exception to be freed.
      doesRethrow(fc.getTarget())
    )
    or
    // Calling operator delete on the exception will free it.
    exists(DeleteExpr d | d = expr | d.getExpr() = exception.getAnAccess())
  )
}

from CatchBlock cb
where
  cb.getParameter().getType().getUnderlyingType() instanceof PointerType and
  not exists(Expr e | e.getEnclosingBlock().getParent*() = cb |
    deletesException(e, cb.getParameter())
  )
select cb, "This catch block does not free the caught exception, thereby leaking memory."
