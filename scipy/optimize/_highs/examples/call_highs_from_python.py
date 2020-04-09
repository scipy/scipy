file = "../src/interfaces/highs_lp_solver.py"
exec(compile(open(file).read(), file, 'exec'))

cc = (1.0, -2.0)
cl = (0.0, 0.0)
cu = (10.0, 10.0)
ru = (2.0, 1.0)
rl = (0.0, 0.0)
astart = (0, 2, 4)
aindex = (0, 1, 0, 1)
avalue = (1.0, 2.0, 1.0, 3.0)
rc, cv, cd, rv, rd , cbs, rbs = Highs_call(cc, cl, cu, rl, ru, astart, aindex, avalue)

print(rc, cv, cd, rv, rd, cbs, rbs)
