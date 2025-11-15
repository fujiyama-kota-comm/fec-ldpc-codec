/**
 * @file ldpc_matrix.h
 * @brief LDPC parity-check (H) and generator (G) matrix utilities.
 *
 * This library provides:
 *  - Regular LDPC H-matrix generation using Gallager construction
 *  - Systematic G-matrix generation from H using Gaussian elimination (GF(2))
 *  - 4-cycle counting utility for LDPC structure analysis
 *
 * All operations are in GF(2), i.e., arithmetic is XOR-based.
 */

#ifndef LDPC_MATRIX_H
#define LDPC_MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================ */
/* 1. LDPC Parity-Check Matrix Generation (Gallager Method)         */
/* ================================================================ */
/**
 * @brief Generate a regular (wc, wr) LDPC parity-check matrix H.
 *
 * @param H  Parity-check matrix (allocated externally): size M × N,
 *           where M = N * wc / wr
 * @param N  Codeword length
 * @param wc Column weight (number of 1s per column)
 * @param wr Row weight (number of 1s per row)
 *
 * Output:
 *   H is filled with 0/1 values using Gallager construction:
 *     - The first block is deterministic
 *     - Remaining (wc-1) blocks are random column permutations
 */
void generate_Hmatrix(int **H, int N, int wc, int wr);

/* ================================================================ */
/* 2. Generator Matrix Construction                                 */
/* ================================================================ */
/**
 * @brief Construct a systematic generator matrix G from an LDPC H matrix.
 *
 * @param H  Parity-check matrix (modified in-place if column swaps occur)
 * @param G  Output generator matrix (allocated externally): size K × N,
 *           where K = N - M and M = N * wc / wr
 * @param N  Codeword length
 * @param wc Column weight
 * @param wr Row weight
 *
 * This function performs the following:
 *   1. Construct extended matrix X = [H^T | I]
 *   2. Perform Gaussian elimination (GF(2))
 *   3. Extract systematic G from bottom K rows of the right block
 *
 * Note:
 *   - Column operations that affect X also update H (left block)
 *   - Produces systematic generator matrix in form G = [P | I]
 */
void generate_Gmatrix(int **H, int **G, int N, int wc, int wr);

/* ================================================================ */
/* 3. Structural Analysis: Count 4-cycles                           */
/* ================================================================ */
/**
 * @brief Count the number of 4-cycles in LDPC parity-check matrix H.
 *
 * @param H  Parity-check matrix (size M × N)
 * @param N  Codeword length
 * @param wc Column weight
 * @param wr Row weight
 *
 * @return Number of 4-cycles detected in H.
 *
 * 4-cycle definition:
 *   Two variable nodes share 2 or more check nodes.
 *
 * Note:
 *   For analysis only. Does not modify H.
 */
int count_floop(int **H, int N, int wc, int wr);

#ifdef __cplusplus
}
#endif

#endif /* LDPC_MATRIX_H */
