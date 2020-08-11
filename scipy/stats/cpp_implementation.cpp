#include <iostream>
#include "cpp_declaration.h"

int hello_world(double *x, int n) {
    for (int i = 0; i < n; i++){
        x[i] = double(i);
    }
    return 0;
}
