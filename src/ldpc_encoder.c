#include "ldpc_encoder.h"

/* =====================================================================
 * LDPC Encoding (Generator-matrix multiplication over GF(2))
 *
 * ecc[i] = Σ_j ( inf[j] · G[j][i] ) mod 2
 *
 * i.e., XOR accumulation:
 *      acc ^= (inf[j] & G[j][i]);
 * ===================================================================== */
void ldpc_encode(int *ecc, const int *inf, int **G, int N, int K) {
  for (int i = 0; i < N; i++) {
    int acc = 0;

    for (int j = 0; j < K; j++) {
      acc ^= (inf[j] & G[j][i]);
    }

    ecc[i] = acc;
  }
}
