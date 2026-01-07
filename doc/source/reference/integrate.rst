.. automodule:: scipy.integrate
   :no-members:
   :no-inherited-members:
   :no-special-members:


Concurrency Notes
-----------------

For guidance on thread-safety of integration solvers, see :ref:`scipy_thread_safety`.

**Thread-safe solvers:**
  - :func:`scipy.integrate.ode` with ``method='lsoda'`` (GIL release + :std:ref:`term-thread-local`)
  - :func:`scipy.integrate.solve_ivp` with ``method='LSODA'`` (GIL release + :std:ref:`term-thread-local`)

**NOT thread-safe:**
  - :func:`scipy.integrate.ode` with ``method='vode'`` - should be :std:ref:`term-externally-synchronized`; use separate
    instance per thread or external lock
