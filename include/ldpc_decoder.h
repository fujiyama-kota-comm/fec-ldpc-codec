#ifndef LDPC_DECODER_H
#define LDPC_DECODER_H

/* =====================================================================
 * LDPC Decoder (Sum-Product Algorithm; SPA)
 *
 * Performs iterative message-passing decoding over Tanner graph:
 *   - Inputs:
 *       LLR[N] : channel log-likelihood ratios
 *       H[M×N] : parity-check matrix
 *       M      : number of parity equations (M = N - K)
 *       N      : code length
 *       K      : information length
 *       max_iter : maximum iteration count
 *
 *   - Outputs:
 *       ecc[N] : decoded codeword bits
 *       inf[K] : decoded information bits (systematic part)
 * ===================================================================== */

#ifdef __cplusplus
extern "C" {
#endif

/* Sum-Product LDPC Decoder (SPA) */
void ldpc_decode_spa(double *LLR, int *ecc, int *inf, int **H, int M, int N,
                     int K, int max_iter);

/* =====================================================================
 * Compute bit-wise LLR from symbol-wise likelihoods
 *
 * pyx[E][N] : likelihood p(y | x = symbol_k) for each symbol index k
 * E         : modulation order (e.g., 2 for BPSK, 4 for QPSK, etc.)
 * N         : number of symbols (each symbol carries log2(E) bits)
 * LLR       : output bit LLR array of length N * log2(E)
 *
 * For bit index b of symbol i:
 *   LLR[b + i*logE] = log( Σ_{k:bit=1} pyx[k][i] / Σ_{k:bit=0} pyx[k][i] )
 * ===================================================================== */
void compute_llr_from_pyx(double **pyx, int E, int N, double *LLR);

#ifdef __cplusplus
}
#endif

#endif /* LDPC_DECODER_H */
