#ifndef SCIPY_HIGHS_HCONFIG_H
#define SCIPY_HIGHS_HCONFIG_H


#ifdef __READERLP_READER_HPP__
// SciPy strdup impl for HiGHS reader.cpp
// inspired by https://stackoverflow.com/a/40766163/17117867

#include <cstring>
#include <cstdlib>

char* strdup (const char* s) {
    std::size_t slen = std::strlen(s);
    char* result = (char*)std::malloc(slen + 1);
    if (NULL == result) {
        return NULL;
    }
    std::memcpy(result, s, slen+1);
    return result;
}
#endif  // __READERLP_READER_HPP__

#endif  // SCIPY_HIGHS_HCONFIG_H
