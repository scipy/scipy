include("../src/interfaces/highs_julia_api.jl")

cc = [1.0, -2.0]
cl = [0.0, 0.0]
cu = [10.0, 10.0]
ru = [2.0, 1.0]
rl = [0.0, 0.0]
astart = [0, 2, 4]
aindex = [0, 1, 0, 1]
avalue = [1.0, 2.0, 1.0, 3.0]

model = HighsModel(cc, cl, cu, rl, ru, astart, aindex, avalue)

# solve model in one go
status, solution, basis, modelstatus = Highs_call(model)

# create highs object and work with it
highs = Highs_create()
Highs_loadModel(highs, model)
Highs_run(highs)
Highs_getSolution(highs)
Highs_getBasis(highs)
Highs_getObjectiveValue(highs)
Highs_getIterationCount(highs)
Highs_destroy(highs)

# read model from file
highs = Highs_create()
Highs_readFromFile(highs, "../check/instances/adlittle.mps")
Highs_run(highs)
Highs_getSolution(highs)
Highs_getBasis(highs)
Highs_getObjectiveValue(highs)
Highs_getIterationCount(highs)
Highs_destroy(highs)

# create model from scratch
highs = Highs_create()
Highs_addCols(highs, cc, cl, cu, [], [], [])
Highs_addRows(highs, rl, ru, astart, aindex, avalue)
Highs_run(highs)
Highs_getSolution(highs)
Highs_getBasis(highs)
Highs_getObjectiveValue(highs)
Highs_getIterationCount(highs)
Highs_destroy(highs)
