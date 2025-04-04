# New C API

The solver options and problem definition options moved to a new structs
that must be initialized before use and the results moved also:
```
prima_problem problem;
prima_init_problem(&problem, n);
problem.x0 = x0;
prima_options options;  
prima_init_options(&options);  
options.iprint = PRIMA_MSG_EXIT;  
options.rhoend= 1e-3;  
options.maxfun = 200*n;  
prima_result result;
const int rc = prima_bobyqa(&fun, &problem, &options, &result);
// print results, or copy them, or use them in computations
prima_free_result(&result);
```
Result struct should be freed after use with its dedicated free function. Make sure to
either use the results or copy them into your own structure before freeing.
