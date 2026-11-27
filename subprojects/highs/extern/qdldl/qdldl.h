/*
 * Elimination tree and up-looking LDL^T factorisation kernels, adapted from
 * QDLDL (https://github.com/osqp/qdldl).
 *
 * Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, The OSQP
 * developers
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Modified for HiPO by the HiGHS developers:
 * - the factorisation computes one row of L per call, so that the caller can
 *   regularise the pivot and check the time limit between rows;
 * - D^-1 is stored on the diagonal of L;
 * - factor_row counts the entries of L it scans;
 * - etree has no input checks, does not return the total count, and its
 *   column counts include the diagonal;
 * - HighsInt indices and int64_t pointers for L.
 */
#ifndef HIGHS_QDLDL_H
#define HIGHS_QDLDL_H

#include "HighsQdldl.h"

#endif
