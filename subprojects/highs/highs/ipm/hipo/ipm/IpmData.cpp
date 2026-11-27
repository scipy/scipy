#include "IpmData.h"

namespace hipo {

void IpmData::append() { record.push_back({}); }
IpmIterData& IpmData::back() { return record.back(); }

}  // namespace hipo
