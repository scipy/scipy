// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "info.h"
#include <cstring>              // for memset
#include <map>
#include <string>
#include "control.h"

ipx_info::ipx_info() {
    memset(this, 0, sizeof(ipx_info));
}

namespace ipx {

template <typename T>
static void dump(std::ostream& os, const char* name, T value) {
    os << Textline(std::string("info.")+name) << value << '\n';
}

std::ostream& operator<<(std::ostream& os, const Info& info) {
    dump(os, "status", info.status);
    dump(os, "status_ipm", info.status_ipm);
    dump(os, "status_crossover", info.status_crossover);
    dump(os, "errflag", info.errflag);

    dump(os, "num_var", info.num_var);
    dump(os, "num_constr", info.num_constr);
    dump(os, "num_entries", info.num_entries);

    dump(os, "num_rows_solver", info.num_rows_solver);
    dump(os, "num_cols_solver", info.num_cols_solver);
    dump(os, "num_entries_solver", info.num_entries_solver);

    dump(os, "dualized", info.dualized);
    dump(os, "dense_cols", info.dense_cols);

    dump(os, "dependent_rows", info.dependent_rows);
    dump(os, "dependent_cols", info.dependent_cols);
    dump(os, "rows_inconsistent", info.rows_inconsistent);
    dump(os, "cols_inconsistent", info.cols_inconsistent);
    dump(os, "primal_dropped", info.primal_dropped);
    dump(os, "dual_dropped", info.dual_dropped);

    dump(os, "abs_presidual", sci2(info.abs_presidual));
    dump(os, "abs_dresidual", sci2(info.abs_dresidual));
    dump(os, "rel_presidual", sci2(info.rel_presidual));
    dump(os, "rel_dresidual", sci2(info.rel_dresidual));
    dump(os, "pobjval", sci8(info.pobjval));
    dump(os, "dobjval", sci8(info.dobjval));
    dump(os, "rel_objgap", sci2(info.rel_objgap));
    dump(os, "complementarity", sci2(info.complementarity));
    dump(os, "normx", sci2(info.normx));
    dump(os, "normy", sci2(info.normy));
    dump(os, "normz", sci2(info.normz));

    dump(os, "objval", sci8(info.objval));
    dump(os, "primal_infeas", sci2(info.primal_infeas));
    dump(os, "dual_infeas", sci2(info.dual_infeas));

    dump(os, "iter", info.iter);
    dump(os, "kktiter1", info.kktiter1);
    dump(os, "kktiter2", info.kktiter2);
    dump(os, "basis_repairs", info.basis_repairs);
    dump(os, "updates_start", info.updates_start);
    dump(os, "updates_ipm", info.updates_ipm);
    dump(os, "updates_crossover", info.updates_crossover);
    dump(os, "pushes_crossover", info.pushes_crossover);

    dump(os, "time_total", fix2(info.time_total));
    dump(os, "time_ipm1", fix2(info.time_ipm1));
    dump(os, "time_ipm2", fix2(info.time_ipm2));
    dump(os, "time_starting_basis", fix2(info.time_starting_basis));
    dump(os, "time_crossover", fix2(info.time_crossover));

    dump(os, "time_kkt_factorize", fix2(info.time_kkt_factorize));
    dump(os, "time_kkt_solve", fix2(info.time_kkt_solve));
    dump(os, "time_maxvol", fix2(info.time_maxvol));
    dump(os, "time_cr1", fix2(info.time_cr1));
    dump(os, "time_cr1_AAt", fix2(info.time_cr1_AAt));
    dump(os, "time_cr1_pre", fix2(info.time_cr1_pre));
    dump(os, "time_cr2", fix2(info.time_cr2));
    dump(os, "time_cr2_NNt", fix2(info.time_cr2_NNt));
    dump(os, "time_cr2_B", fix2(info.time_cr2_B));
    dump(os, "time_cr2_Bt", fix2(info.time_cr2_Bt));

    dump(os, "ftran_sparse", fix2(info.ftran_sparse));
    dump(os, "btran_sparse", fix2(info.btran_sparse));
    dump(os, "time_ftran", fix2(info.time_ftran));
    dump(os, "time_btran", fix2(info.time_btran));
    dump(os, "time_lu_invert", fix2(info.time_lu_invert));
    dump(os, "time_lu_update", fix2(info.time_lu_update));
    dump(os, "mean_fill", fix2(info.mean_fill));
    dump(os, "max_fill", fix2(info.max_fill));
    dump(os, "time_symb_invert", fix2(info.time_symb_invert));

    dump(os, "maxvol_updates", info.maxvol_updates);
    dump(os, "maxvol_skipped", info.maxvol_skipped);
    dump(os, "maxvol_passes", info.maxvol_passes);
    dump(os, "tbl_nnz", info.tbl_nnz);
    dump(os, "tbl_max", sci2(info.tbl_max));
    dump(os, "frobnorm_squared", sci2(info.frobnorm_squared));
    dump(os, "lambdamax", sci2(info.lambdamax));
    dump(os, "volume_increase", sci2(info.volume_increase));

    return os;
}

std::string StatusString(Int status) {
    const std::map<int,std::string> status_name = {
        {IPX_STATUS_not_run, "not run"},
        {IPX_STATUS_solved, "solved"},
        {IPX_STATUS_stopped, "stopped"},
        {IPX_STATUS_invalid_input, "invalid input"},
        {IPX_STATUS_out_of_memory, "out of memory"},
        {IPX_STATUS_internal_error, "internal error"},
        {IPX_STATUS_optimal, "optimal"},
        {IPX_STATUS_imprecise, "imprecise"},
        {IPX_STATUS_primal_infeas, "primal infeas"},
        {IPX_STATUS_dual_infeas, "dual infeas"},
        {IPX_STATUS_time_limit, "time limit"},
        {IPX_STATUS_iter_limit, "iter limit"},
        {IPX_STATUS_no_progress, "no progress"},
        {IPX_STATUS_failed, "failed"},
        {IPX_STATUS_debug, "debug"}
    };
    auto entry = status_name.find(status);
    if (entry != status_name.end())
        return entry->second;
    return "unknown";
}

}  // namespace ipx
