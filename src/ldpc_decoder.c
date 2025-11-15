#include "ldpc_decoder.h"
#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------
 * Helper: sign(x)
 * ------------------------------------------------------------ */
static inline double sign_val(double x) { return (x >= 0.0) ? 1.0 : -1.0; }

/* ------------------------------------------------------------
 * Helper: spf(x) = log((e^x + 1)/(e^x - 1))
 * With safe limits to avoid overflow / underflow
 * ------------------------------------------------------------ */
static inline double spf(double x) {
  if (x < 1e-7)
    x = 1e-7;
  else if (x > 30.0)
    x = 30.0;

  return log((exp(x) + 1.0) / (exp(x) - 1.0));
}

/* ============================================================
 * Sum-Product LDPC Decoder (SPA)
 * ============================================================ */
void ldpc_decode_spa(double *LLR, int *ecc, int *inf, int **H, int M, int N,
                     int K, int max_iter) {
  int i, j, k, iter;

  /* ------------------------------------
   * Build check node and variable node lists
   * ------------------------------------ */
  int **check_node = (int **)malloc(M * sizeof(int *));
  for (i = 0; i < M; i++) {
    check_node[i] = (int *)malloc(N * sizeof(int));
  }

  int *deg_c = (int *)calloc(M, sizeof(int));
  int *deg_v = (int *)calloc(N, sizeof(int));

  /* Count degrees of check nodes and variable nodes */
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++) {
      if (H[i][j])
        deg_c[i]++;
    }
  }
  for (j = 0; j < N; j++) {
    for (i = 0; i < M; i++) {
      if (H[i][j])
        deg_v[j]++;
    }
  }

  /* Variable-node adjacency list */
  int **variable_node = (int **)malloc(N * sizeof(int *));
  for (j = 0; j < N; j++) {
    variable_node[j] = (int *)malloc(deg_v[j] * sizeof(int));
  }

  /* build check_node list */
  for (i = 0; i < M; i++) {
    int idx = 0;
    for (j = 0; j < N; j++) {
      if (H[i][j]) {
        check_node[i][idx++] = j;
      }
    }
  }

  /* build variable_node list */
  for (j = 0; j < N; j++) {
    int idx = 0;
    for (i = 0; i < M; i++) {
      if (H[i][j]) {
        variable_node[j][idx++] = i;
      }
    }
  }

  /* ------------------------------------
   * Allocate messages u[i][j] and v[i][j]
   * ------------------------------------ */
  double **u = (double **)malloc(M * sizeof(double *));
  double **v = (double **)malloc(M * sizeof(double *));
  for (i = 0; i < M; i++) {
    u[i] = (double *)calloc(N, sizeof(double));
    v[i] = (double *)calloc(N, sizeof(double));
  }

  /* ============================================================
   * Iterative SPA
   * ============================================================ */
  for (iter = 0; iter < max_iter; iter++) {

    /* ---------------------------
     * Check node update
     * --------------------------- */
    for (i = 0; i < M; i++) {
      for (k = 0; k < deg_c[i]; k++) {

        int j_idx = check_node[i][k];

        double prod_sign = 1.0;
        double sum_spf_val = 0.0;

        for (j = 0; j < deg_c[i]; j++) {
          int var = check_node[i][j];
          if (j != k) {
            double x = LLR[var] + u[i][var];
            prod_sign *= sign_val(x);
            sum_spf_val += spf(fabs(x));
          }
        }

        v[i][j_idx] = prod_sign * spf(sum_spf_val);
      }
    }

    /* ---------------------------
     * Variable node update
     * --------------------------- */
    for (j = 0; j < N; j++) {
      for (k = 0; k < deg_v[j]; k++) {

        int i_idx = variable_node[j][k];
        double sum_v = 0.0;

        for (i = 0; i < deg_v[j]; i++) {
          int cnode = variable_node[j][i];
          if (i != k) {
            sum_v += v[cnode][j];
          }
        }
        u[i_idx][j] = sum_v;
      }
    }

    /* ---------------------------
     * Compute tentative decision
     * --------------------------- */
    for (j = 0; j < N; j++) {
      double sum = LLR[j];
      for (i = 0; i < deg_v[j]; i++) {
        sum += v[variable_node[j][i]][j];
      }
      ecc[j] = (sum >= 0.0) ? 1 : 0;
    }

    /* ---------------------------
     * Parity check H * ecc^T == 0 ?
     * --------------------------- */
    int parity_ok = 1;
    for (i = 0; i < M; i++) {
      int parity = 0;
      for (k = 0; k < deg_c[i]; k++) {
        parity ^= ecc[check_node[i][k]];
      }
      if (parity != 0) {
        parity_ok = 0;
        break;
      }
    }

    if (parity_ok) {
      break; /* successful decoding */
    }
  }

  /* ------------------------------------
   * Extract information bits (systematic part)
   *   Assumes codeword = [parity | info]
   * ------------------------------------ */
  for (i = 0; i < K; i++) {
    inf[i] = ecc[i + (N - K)];
  }

  /* ------------------------------------
   * Free memory
   * ------------------------------------ */
  for (i = 0; i < M; i++) {
    free(check_node[i]);
  }
  free(check_node);

  for (j = 0; j < N; j++) {
    free(variable_node[j]);
  }
  free(variable_node);

  free(deg_c);
  free(deg_v);

  for (i = 0; i < M; i++) {
    free(u[i]);
    free(v[i]);
  }
  free(u);
  free(v);
}

/* ============================================================
 * Compute LLR from pyx[E][N]
 * ============================================================ */
void compute_llr_from_pyx(double **pyx, int E, int N, double *LLR) {
  int i, k, b;
  int logE = (int)log2((double)E);

  /* Convert symbol index → bits */
  int **symbol_bits = (int **)malloc(E * sizeof(int *));
  for (k = 0; k < E; k++) {
    symbol_bits[k] = (int *)malloc(logE * sizeof(int));
    for (b = 0; b < logE; b++) {
      symbol_bits[k][b] = (k >> b) & 1;
    }
  }

  /* Compute bit-wise LLR */
  for (i = 0; i < N; i++) {
    for (b = 0; b < logE; b++) {

      double p1 = 0.0;
      double p0 = 0.0;

      for (k = 0; k < E; k++) {
        if (symbol_bits[k][b])
          p1 += pyx[k][i];
        else
          p0 += pyx[k][i];
      }

      LLR[b + i * logE] = log(p1 / p0);
    }
  }

  for (k = 0; k < E; k++) {
    free(symbol_bits[k]);
  }
  free(symbol_bits);
}
