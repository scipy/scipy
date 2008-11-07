/**
 * common.h
 *
 * Author: Damian Eads
 * Date:   September 22, 2007 (moved into new file on June 8, 2008)
 *
 * Copyright (c) 2007, 2008, Damian Eads. All rights reserved.
 * Adapted for incorporation into Scipy, April 9, 2008.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   - Redistributions of source code must retain the above
 *     copyright notice, this list of conditions and the
 *     following disclaimer.
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   - Neither the name of the author nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _CLUSTER_COMMON_H
#define _CLUSTER_COMMON_H

#define CPY_MAX(_x, _y) ((_x > _y) ? (_x) : (_y))
#define CPY_MIN(_x, _y) ((_x < _y) ? (_x) : (_y))

#define NCHOOSE2(_n) ((_n)*(_n-1)/2)

#define CPY_BITS_PER_CHAR (sizeof(unsigned char) * 8)
#define CPY_FLAG_ARRAY_SIZE_BYTES(num_bits) (CPY_CEIL_DIV((num_bits), \
                                                          CPY_BITS_PER_CHAR))
#define CPY_GET_BIT(_xx, i) (((_xx)[(i) / CPY_BITS_PER_CHAR] >> \
                             ((CPY_BITS_PER_CHAR-1) - \
                              ((i) % CPY_BITS_PER_CHAR))) & 0x1)
#define CPY_SET_BIT(_xx, i) ((_xx)[(i) / CPY_BITS_PER_CHAR] |= \
                              ((0x1) << ((CPY_BITS_PER_CHAR-1) \
                                         -((i) % CPY_BITS_PER_CHAR))))
#define CPY_CLEAR_BIT(_xx, i) ((_xx)[(i) / CPY_BITS_PER_CHAR] &= \
                              ~((0x1) << ((CPY_BITS_PER_CHAR-1) \
                                         -((i) % CPY_BITS_PER_CHAR))))

#ifndef CPY_CEIL_DIV
#define CPY_CEIL_DIV(x, y) ((((double)x)/(double)y) == \
                            ((double)((x)/(y))) ? ((x)/(y)) : ((x)/(y) + 1))
#endif

#ifdef CPY_DEBUG
#define CPY_DEBUG_MSG(...) fprintf(stderr, __VA_ARGS__)
#else
#define CPY_DEBUG_MSG(...)
#endif

#endif
