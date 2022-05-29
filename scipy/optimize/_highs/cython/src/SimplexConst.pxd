# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool

cdef extern from "SimplexConst.h" nogil:

    cdef enum SimplexAlgorithm:
        PRIMAL "SimplexAlgorithm::kPrimal" = 0
        DUAL "SimplexAlgorithm::kDual"

    cdef enum SimplexStrategy:
        SIMPLEX_STRATEGY_MIN "SimplexStrategy::kSimplexStrategyMin" = 0
        SIMPLEX_STRATEGY_CHOOSE "SimplexStrategy::kSimplexStrategyChoose" = SIMPLEX_STRATEGY_MIN
        SIMPLEX_STRATEGY_DUAL "SimplexStrategy::kSimplexStrategyDual"
        SIMPLEX_STRATEGY_DUAL_PLAIN "SimplexStrategy::kSimplexStrategyDualPlain" = SIMPLEX_STRATEGY_DUAL
        SIMPLEX_STRATEGY_DUAL_TASKS "SimplexStrategy::kSimplexStrategyDualTasks"
        SIMPLEX_STRATEGY_DUAL_MULTI "SimplexStrategy::kSimplexStrategyDualMulti"
        SIMPLEX_STRATEGY_PRIMAL "SimplexStrategy::kSimplexStrategyPrimal"
        SIMPLEX_STRATEGY_MAX "SimplexStrategy::kSimplexStrategyMax" = SIMPLEX_STRATEGY_PRIMAL
        SIMPLEX_STRATEGY_NUM "SimplexStrategy::kSimplexStrategyNum"

    cdef enum SimplexCrashStrategy:
        SIMPLEX_CRASH_STRATEGY_MIN "SimplexCrashStrategy::kSimplexCrashStrategyMin" = 0
        SIMPLEX_CRASH_STRATEGY_OFF "SimplexCrashStrategy::kSimplexCrashStrategyOff" = SIMPLEX_CRASH_STRATEGY_MIN
        SIMPLEX_CRASH_STRATEGY_LTSSF_K "SimplexCrashStrategy::kSimplexCrashStrategyLtssfK"
        SIMPLEX_CRASH_STRATEGY_LTSSF "SimplexCrashStrategy::kSimplexCrashStrategyLtssf" = SIMPLEX_CRASH_STRATEGY_LTSSF_K
        SIMPLEX_CRASH_STRATEGY_BIXBY "SimplexCrashStrategy::kSimplexCrashStrategyBixby"
        SIMPLEX_CRASH_STRATEGY_LTSSF_PRI "SimplexCrashStrategy::kSimplexCrashStrategyLtssfPri"
        SIMPLEX_CRASH_STRATEGY_LTSF_K "SimplexCrashStrategy::kSimplexCrashStrategyLtsfK"
        SIMPLEX_CRASH_STRATEGY_LTSF_PRI "SimplexCrashStrategy::kSimplexCrashStrategyLtsfPri"
        SIMPLEX_CRASH_STRATEGY_LTSF "SimplexCrashStrategy::kSimplexCrashStrategyLtsf"
        SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS "SimplexCrashStrategy::kSimplexCrashStrategyBixbyNoNonzeroColCosts"
        SIMPLEX_CRASH_STRATEGY_BASIC "SimplexCrashStrategy::kSimplexCrashStrategyBasic"
        SIMPLEX_CRASH_STRATEGY_TEST_SING "SimplexCrashStrategy::kSimplexCrashStrategyTestSing"
        SIMPLEX_CRASH_STRATEGY_MAX "SimplexCrashStrategy::kSimplexCrashStrategyMax" = SIMPLEX_CRASH_STRATEGY_TEST_SING

    cdef enum SimplexDualEdgeWeightStrategy:
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategyMin" = -1
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_CHOOSE "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategyChoose" = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategyDantzig"
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategyDevex"
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategySteepestEdge"
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategySteepestEdgeUnitInitial"
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MAX "SimplexDualEdgeWeightStrategy::kSimplexDualEdgeWeightStrategyMax" = SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL

    cdef enum SimplexPrimalEdgeWeightStrategy:
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN "SimplexPrimalEdgeWeightStrategy::kSimplexPrimalEdgeWeightStrategyMin" = -1
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_CHOOSE "SimplexPrimalEdgeWeightStrategy::kSimplexPrimalEdgeWeightStrategyChoose" = SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG "SimplexPrimalEdgeWeightStrategy::kSimplexPrimalEdgeWeightStrategyDantzig"
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX "SimplexPrimalEdgeWeightStrategy::kSimplexPrimalEdgeWeightStrategyDevex"
        SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MAX "SimplexPrimalEdgeWeightStrategy::kSimplexPrimalEdgeWeightStrategyMax" = SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX

    cdef enum SimplexPriceStrategy:
        SIMPLEX_PRICE_STRATEGY_MIN = 0
        SIMPLEX_PRICE_STRATEGY_COL = SIMPLEX_PRICE_STRATEGY_MIN
        SIMPLEX_PRICE_STRATEGY_ROW
        SIMPLEX_PRICE_STRATEGY_ROW_SWITCH
        SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH
        SIMPLEX_PRICE_STRATEGY_MAX = SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH

    cdef enum SimplexDualChuzcStrategy:
        SIMPLEX_DUAL_CHUZC_STRATEGY_MIN = 0
        SIMPLEX_DUAL_CHUZC_STRATEGY_CHOOSE = SIMPLEX_DUAL_CHUZC_STRATEGY_MIN
        SIMPLEX_DUAL_CHUZC_STRATEGY_QUAD
        SIMPLEX_DUAL_CHUZC_STRATEGY_HEAP
        SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH
        SIMPLEX_DUAL_CHUZC_STRATEGY_MAX = SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH

    cdef enum InvertHint:
        INVERT_HINT_NO = 0
        INVERT_HINT_UPDATE_LIMIT_REACHED
        INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT
        INVERT_HINT_POSSIBLY_OPTIMAL
        INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED
        INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED
        INVERT_HINT_POSSIBLY_SINGULAR_BASIS
        INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX
        INVERT_HINT_CHOOSE_COLUMN_FAIL
        INVERT_HINT_Count

    cdef enum DualEdgeWeightMode:
        DANTZIG "DualEdgeWeightMode::DANTZIG" = 0
        DEVEX "DualEdgeWeightMode::DEVEX"
        STEEPEST_EDGE "DualEdgeWeightMode::STEEPEST_EDGE"
        Count "DualEdgeWeightMode::Count"

    cdef enum PriceMode:
        ROW "PriceMode::ROW" = 0
        COL "PriceMode::COL"

    const int PARALLEL_THREADS_DEFAULT
    const int DUAL_TASKS_MIN_THREADS
    const int DUAL_MULTI_MIN_THREADS

    const bool invert_if_row_out_negative

    const int NONBASIC_FLAG_TRUE
    const int NONBASIC_FLAG_FALSE

    const int NONBASIC_MOVE_UP
    const int NONBASIC_MOVE_DN
    const int NONBASIC_MOVE_ZE
