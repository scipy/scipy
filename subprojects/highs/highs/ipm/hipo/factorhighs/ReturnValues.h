#ifndef FACTORHIGHS_RETURN_VALUES_H
#define FACTORHIGHS_RETURN_VALUES_H

namespace hipo {

enum RetValue {
  kRetOk,
  kRetInvalidInput,
  kRetOutOfMemory,
  kRetInvalidPivot,
  kRetOrderingError,
  kRetIntOverflow,
  kRetInvalidPointer,
  kRetGeneric
};

}

#endif