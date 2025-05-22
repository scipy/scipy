#ifndef __DIA_H__
#define __DIA_H__

#include <algorithm>
#include <vector>

using std::min, std::max;

template <class I, class T>
T min_el(const T vec[], I len) {
  return *std::min_element(vec, vec + len);
}

template <class I, class T>
T max_el(const T vec[], I len) {
  return *std::max_element(vec, vec + len);
}


/*
 * Compute DIA matrix output = A * B for DIA matrices A, B
 *
 *
 * Input Arguments:
 *   I  A_rows               - number of rows in A
 *   I  A_cols               - number of columns in A
 *   I  A_diags              - number of diagonals in A
 *   I  A_L                  - length of each diagonal in A
 *   I  A_offsets[A_diags]   - diagonal offsets in A
 *   T  A_data[A_diags,A_L]  - diagonals data of A
 *   I  B_cols               - number of columns in B
 *   I  B_diags              - number of diagonals in B
 *   I  B_L                  - length of each diagonal in B
 *   I  B_offsets[B_diags]   - diagonal offsets in B
 *   T  B_data[B_diags,B_L]  - diagonals data of B
 *
 * Output Arguments:
 *   V  offsets              - diagonal offsets
 *   W  data                 - diagonals data
 *
 * Note:
 *   Number of diagonals in output (diags) is the length of output offsets;
 *   length of each diagonal (L) equals B_L; data dimensions are [diags,L]
 *   Negative offsets correspond to lower diagonals
 *   Positive offsets correspond to upper diagonals
 *
 */
template <class I, class T>
void dia_matmat(const I A_rows,
                const I A_cols,
                const I A_diags,
                const I A_L,
                const I A_offsets[],
                const T A_data[],
                const I B_cols,
                const I B_diags,
                const I B_L,
                const I B_offsets[],
                const T B_data[],
                std::vector<I>* offsets,
                std::vector<T>* data)
{
    const I B_rows = A_cols,
            rows = A_rows, cols = B_cols,
            L = min(cols, B_L);

    // range of combinations of input offsets
    const I min_map_ofs = min_el(A_offsets, A_diags) + min_el(B_offsets, B_diags),
            max_map_ofs = max_el(A_offsets, A_diags) + max_el(B_offsets, B_diags);
    // limits for output offsets
    const I min_ofs = max(min_map_ofs, 1 - rows),
            max_ofs = min(max_map_ofs, L - 1);
    // mark output offsets
    std::vector<I> buf(max_map_ofs - min_map_ofs + 1);
    // (alias to buf indexable from [min_map_ofs] to [max_map_ofs];
    //  only [min_ofs] to [max_ofs] will be used later)
    I* ofs_map = buf.data() - min_map_ofs;
    for (I i = 0; i < A_diags; ++i)
        for (I j = 0; j < B_diags; ++j)
            ofs_map[A_offsets[i] + B_offsets[j]] = I(true);
    // enumerate marks
    offsets->resize(max_ofs - min_ofs + 1);
    I N_ofs = 0;
    for (I ofs = min_ofs; ofs <= max_ofs; ++ofs)
        if (ofs_map[ofs]) {
            (*offsets)[N_ofs] = ofs;
            ofs_map[ofs] = N_ofs;
            ++N_ofs;
        }
    offsets->resize(N_ofs);

    // allocate output diagonals, filled with zeros
    data->resize(N_ofs * L);
    // loop over diagonals in B
    for (I B_i = 0; B_i < B_diags; ++B_i) {
        const I B_ofs = B_offsets[B_i];
        const T* B_diag_r = B_data + npy_intp(B_L) * B_i + B_ofs; // row-indexed
        // span of data within current (row-indexed) B diagonal
        const I B_j_beg = max(-B_ofs, I(0)),
                B_j_end = min(min(B_cols, B_L) - B_ofs, B_rows);
        // loop over diagonals in A
        for (I A_i = 0; A_i < A_diags; ++A_i) {
            const I A_ofs = A_offsets[A_i];
            const T* A_diag_c = A_data + npy_intp(A_L) * A_i; // column-indexed
            // output offset and corresponding diagonal
            const I ofs = A_ofs + B_ofs;
            if (ofs < min_ofs or ofs > max_ofs)
                continue;
            T* diag_r = data->data() + npy_intp(L) * ofs_map[ofs] + B_ofs; // row-indexed
            // overlapping span of data within current B and A diagonals
            const I j_beg = max(B_j_beg, A_ofs),
                    j_end = min({B_j_end, min(A_cols, A_L), A_rows + A_ofs});
            // add partial product to output
            for (I j = j_beg; j < j_end; ++j)
                diag_r[j] += A_diag_c[j] * B_diag_r[j];
        }
    }
}


/*
 * Compute Y += A * X for DIA matrix A and dense vectors X, Y
 *
 *
 * Input Arguments:
 *   I  n_row            - number of rows in A
 *   I  n_col            - number of columns in A
 *   I  n_diags          - number of diagonals
 *   I  L                - length of each diagonal
 *   I  offsets[n_diags] - diagonal offsets 
 *   T  diags[n_diags,L] - nonzeros 
 *   T  Xx[n_col]        - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]        - output vector 
 *
 * Note:
 *   Output array Yx must be preallocated
 *   Negative offsets correspond to lower diagonals
 *   Positive offsets correspond to upper diagonals
 *
 */
template <class I, class T>
void dia_matvec(const I n_row,
                const I n_col,
                const I n_diags,
                const I L,
                const I offsets[],
                const T diags[],
                const T Xx[],
                      T Yx[])
{
    for (I i = 0; i < n_diags; i++){
        const I k = offsets[i];  //diagonal offset

        const I i_start = max<I>(0, -k);
        const I j_start = max<I>(0, k);
        const I j_end   = min<I>(min<I>(n_row + k, n_col), L);

        const I N = j_end - j_start;  //number of elements to process

        const T * diag = diags + (npy_intp)i * L + j_start;
        const T * x = Xx + j_start;
              T * y = Yx + i_start;

        for (I n = 0; n < N; n++) {
            y[n] += diag[n] * x[n]; 
        }
    }
}


/*
 * Compute output += A * B for DIA matrix A and dense matrices B, output
 *
 *
 * Input Arguments:
 *   I  A_rows                - number of rows in A
 *   I  A_cols                - number of columns in A
 *   I  A_diags               - number of diagonals in A
 *   I  A_L                   - length of each diagonal in A
 *   I  offsets[A_diags]      - diagonal offsets in A
 *   T  A_data[A_diags,A_L]   - diagonals data of A
 *   I  B_cols                - number of columns in B
 *   T  B_data[A_rows,B_cols] - data of B (in C order)
 *
 * Output Arguments:
 *   T  data[A_rows,B_cols]   - output data (in C order)
 *
 * Note:
 *   Output array data must be preallocated
 *   Number of rows in B must be equal A_cols
 *   Negative offsets correspond to lower diagonals
 *   Positive offsets correspond to upper diagonals
 *
 */
template <class I, class T>
void dia_matvecs(const I A_rows,
                 const I A_cols,
                 const I A_diags,
                 const I A_L,
                 const I A_offsets[],
                 const T A_data[],
                 const I B_cols,
                 const T B_data[],
                       T data[])
{
    const I rows = A_rows, cols = B_cols,
            k_end = min(A_cols, A_L); // for index along A columns and B rows
    // loop over output rows
    for (I i = 0; i < rows; ++i) {
        T* row = data + npy_intp(cols) * i;
        // loop over diagonals in A
        for (I n = 0; n < A_diags; ++n) {
            const I k = i + A_offsets[n];
            if (k < 0 or k >= k_end)
                continue;
            // element at i-th row, k-th column in A
            const T a = (A_data + npy_intp(A_L) * n)[k];
            // k-th row in B
            const T* B_row = B_data + npy_intp(B_cols) * k;
            // loop over columns in current output row
            for (I j = 0; j < cols; ++j)
                row[j] += a * B_row[j];
        }
    }
}


#endif
