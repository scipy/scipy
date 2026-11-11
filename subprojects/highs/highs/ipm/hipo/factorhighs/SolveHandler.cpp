#include "SolveHandler.h"

namespace hipo {

SolveHandler::SolveHandler(const Symbolic& S,
                           const std::vector<std::vector<double>>& sn_columns,
                           DataCollector& data, const FHoptions& options)
    : S_{S}, sn_columns_{sn_columns}, data_{data}, options_{options} {}

}  // namespace hipo