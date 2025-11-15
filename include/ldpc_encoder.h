#ifndef LDPC_ENCODER_H
#define LDPC_ENCODER_H

/* =====================================================================
 * LDPC Encoder (Generator-matrix based)
 *
 * Performs:
 *      c = u × G   over GF(2)
 *
 * Inputs:
 *      inf[K]  : information bits   (0/1)
 *      G[K][N] : generator matrix
 *      N       : codeword length
 *      K       : information length
 *
 * Output:
 *      ecc[N]  : encoded LDPC codeword
 * ===================================================================== */

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------
 * ecc[i] = XOR_j (inf[j] & G[j][i])
 * GF(2) linear encoding
 * --------------------------------------------------------- */
void ldpc_encode(int *ecc, const int *inf, int **G, int N, int K);

#ifdef __cplusplus
}
#endif

#endif /* LDPC_ENCODER_H */
