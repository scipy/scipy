#ifndef RCONT_H
#define RCONT_H

#include <stdint.h>

#include "distributions.h"

typedef int64_t tab_t;

void rcont1_init(tab_t *work, int nc, const tab_t *c);

void rcont1(tab_t *table, int nr, const tab_t *r, int nc, const tab_t *c,
            const tab_t ntot, tab_t *work, bitgen_t *rstate);

void rcont2(tab_t *table, int nr, const tab_t *r, int nc, const tab_t *c,
            const tab_t ntot, bitgen_t *rstate);

#endif