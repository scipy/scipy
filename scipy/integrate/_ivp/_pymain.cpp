// Compile riccaticpp's pybind11 bindings as part of scipy.
//
// meson's sandbox forbids files() references that cross into a
// subproject directory, so we cannot list
// ``subprojects/riccaticpp/src/pymain.cpp`` directly in
// scipy/integrate/_ivp/meson.build. Including it here lets the C++
// preprocessor pull the upstream source verbatim -- no duplicated
// code, no drift, just one line of glue.
#include "../../../subprojects/riccaticpp/src/pymain.cpp"
