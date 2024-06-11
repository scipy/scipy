#pragma once

#if defined(SCIPY_DLL_EXPORTS)
    #define SCIPY_DLL __declspec(dllexport)
#elif defined(SCIPY_DLL_IMPORTS)
    #define SCIPY_DLL __declspec(dllimport)
#else
    #define SCIPY_DLL
#endif
