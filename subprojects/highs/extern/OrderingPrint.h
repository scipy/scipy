#ifndef HIPO_ORDERING_PRINT_H
#define HIPO_ORDERING_PRINT_H

#ifndef NDEBUG
#define HIGHS_ORDERING_PRINT(...) printf(__VA_ARGS__)
#else
#define HIGHS_ORDERING_PRINT(...)
#endif

#endif