#pragma once

// SCIPY_DLL
// inspired by https://github.com/abseil/abseil-cpp/blob/20240116.2/absl/base/config.h#L736-L753
//
// When building sf_error_state as a DLL, this macro expands to `__declspec(dllexport)`
// so we can annotate symbols appropriately as being exported. When used in
// headers consuming a DLL, this macro expands to `__declspec(dllimport)` so
// that consumers know the symbol is defined inside the DLL. In all other cases,
// the macro expands to nothing.
// Note: SCIPY_DLL_{EX,IM}PORTS are set in scipy/special/meson.build
#if defined(SCIPY_DLL_EXPORTS)
    #define SCIPY_DLL __declspec(dllexport)
#elif defined(SCIPY_DLL_IMPORTS)
    #define SCIPY_DLL __declspec(dllimport)
#else
    #define SCIPY_DLL
#endif
