/**
 * @file ldpc_matrix.c
 * @brief LDPC parity-check matrix (H) and generator matrix (G) utilities
 *
 * This file contains:
 *   - Gallager-type regular LDPC H-matrix generation
 *   - Systematic generator matrix G construction (via Gaussian elimination)
 *   - 4-cycle counting utility (analysis of LDPC structure)
 *
 * すべて行列演算は GF(2)（XOR）で行われます。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ldpc_matrix.h"

/* ================================================================
 * 1. Generate regular (w_c, w_r) LDPC parity-check matrix (Gallager)
 * ----------------------------------------------------------------
 * Output: H[M][N] where M = N*(w_c/w_r)
 * H は (w_c × サブブロック数) 個の行ブロックで構成される。
 * 第一ブロックは連続した1、以降は列のランダムパーミュテーションで生成。
 * ================================================================ */
void generate_Hmatrix(int **H, int N, int wc, int wr) {
  int i, j, k;
  int M = (N * wc) / wr;   /* number of parity equations */
  int block_rows = M / wc; /* rows per block (Gallager construction) */

  /* permutation buffer */
  int *perm = (int *)malloc(N * sizeof(int));
  if (!perm) {
    fprintf(stderr, "malloc failed in generate_Hmatrix\n");
    exit(1);
  }

  /* Clear H matrix */
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++)
      H[i][j] = 0;

  /**
   * Step 1: First block (identity-type pattern)
   * Rows: 0 ... block_rows-1
   * Columns: i*wr ... (i+1)*wr-1
   */
  for (i = 0; i < block_rows; i++) {
    for (j = i * wr; j < (i + 1) * wr; j++) {
      H[i][j] = 1;
    }
  }

  /**
   * Step 2: Remaining block rows are permutations of the first block
   */
  for (i = 1; i < wc; i++) {

    /* generate identity permutation */
    for (j = 0; j < N; j++)
      perm[j] = j;

    /* shuffle columns randomly */
    for (j = 0; j < N; j++) {
      int r = rand() % N;
      int tmp = perm[j];
      perm[j] = perm[r];
      perm[r] = tmp;
    }

    /* apply permuted columns over next block */
    for (j = block_rows * i; j < block_rows * (i + 1); j++) {
      for (k = 0; k < N; k++) {
        H[j][k] = H[j - block_rows * i][perm[k]];
      }
    }
  }

  free(perm);
}

/* ================================================================
 * 2. Construct systematic generator matrix G from H
 * ----------------------------------------------------------------
 * Input:  H[M][N]
 * Output: G[K][N] where K = N - M
 *
 * 方法：
 *   [H^T | I_N] を拡張行列としてガウス消去し、
 *   下側 K 行から右側 N 列部分を取り出して G を作る。
 *
 * 注意：元のアルゴリズムは温存（行・列交換も含む）
 * ================================================================ */
void generate_Gmatrix(int **H, int **G, int N, int wc, int wr) {
  int M = (N * wc) / wr;
  int K = N - M;

  int i, j, k, l;

  /* X: transformation matrix of size N × (M+N) */
  int **X = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++)
    X[i] = (int *)malloc((M + N) * sizeof(int));

  int *row_buf = (int *)malloc((M + N) * sizeof(int));
  int *col_buf = (int *)malloc(N * sizeof(int));
  int *Hcol_buf = (int *)malloc(M * sizeof(int));

  if (!X || !row_buf || !col_buf || !Hcol_buf) {
    fprintf(stderr, "malloc failed in generate_Gmatrix\n");
    exit(1);
  }

  /* ------------------------------------------------------------
   * Step 1: Construct extended matrix [H^T | I]
   * ------------------------------------------------------------ */
  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++)
      X[i][j] = H[j][i]; /* Left block = H^T */

    for (j = M; j < M + N; j++)
      X[i][j] = (i == (j - M)) ? 1 : 0; /* Right block = identity */
  }

  /* ------------------------------------------------------------
   * Step 2: Gaussian elimination on the left part of X
   *         (Column swaps here DO NOT affect H)
   * ------------------------------------------------------------ */
  for (j = 0; j < M; j++) {
    if (X[j][j] == 0) {
      /* Find pivot row */
      int pivot_found = 0;
      for (i = j + 1; i < N; i++) {
        if (X[i][j] == 1) {
          memcpy(row_buf, X[i], (M + N) * sizeof(int));
          memcpy(X[i], X[j], (M + N) * sizeof(int));
          memcpy(X[j], row_buf, (M + N) * sizeof(int));
          pivot_found = 1;
          break;
        }
      }
      /* If pivot row not found, swap columns */
      if (!pivot_found) {
        for (k = M + N - 1; k > j; k--) {
          if (X[j][k] == 1) {
            for (l = 0; l < N; l++) {
              col_buf[l] = X[l][k];
              X[l][k] = X[l][j];
              X[l][j] = col_buf[l];
            }
            break;
          }
        }
      }
    }

    /* row elimination */
    for (i = 0; i < N; i++) {
      if (i != j && X[i][j] == 1) {
        for (k = 0; k < (M + N); k++)
          X[i][k] ^= X[j][k]; /* XOR */
      }
    }
  }

  /* ------------------------------------------------------------
   * Step 3: Elimination that affects both X and H (column swaps)
   * ------------------------------------------------------------ */
  for (j = 2 * M; j < M + N; j++) {

    int pivot_row = j - M;

    if (X[pivot_row][j] == 0) {

      /* find pivot row */
      int found = 0;
      for (i = pivot_row + 1; i < N; i++) {
        if (X[i][j] == 1) {
          memcpy(row_buf, X[i], (M + N) * sizeof(int));
          memcpy(X[i], X[pivot_row], (M + N) * sizeof(int));
          memcpy(X[pivot_row], row_buf, (M + N) * sizeof(int));
          found = 1;
          break;
        }
      }

      /* pivot row not found ⇒ column swap (affects H also) */
      if (!found) {
        for (k = (M + N) - 1; k > M - 1; k--) {
          if (X[pivot_row][k] == 1) {

            /* swap columns in X */
            for (l = 0; l < N; l++) {
              col_buf[l] = X[l][k];
              X[l][k] = X[l][j];
              X[l][j] = col_buf[l];
            }

            /* swap corresponding columns in H */
            for (l = 0; l < M; l++) {
              Hcol_buf[l] = H[l][k - M];
              H[l][k - M] = H[l][j - M];
              H[l][j - M] = Hcol_buf[l];
            }

            break;
          }
        }
      }
    }

    /* eliminate other rows */
    for (i = 0; i < N; i++) {
      if (i != pivot_row && X[i][j] == 1) {
        for (k = 0; k < M + N; k++)
          X[i][k] ^= X[pivot_row][k];
      }
    }
  }

  /* ------------------------------------------------------------
   * Step 4: Extract generator matrix G  (K × N block)
   * G = X[M...N-1][M...M+N-1]
   * ------------------------------------------------------------ */
  for (i = M; i < N; i++)
    for (j = M; j < M + N; j++)
      G[i - M][j - M] = X[i][j];

  /* free memory */
  for (i = 0; i < N; i++)
    free(X[i]);
  free(X);
  free(row_buf);
  free(col_buf);
  free(Hcol_buf);
}

/* ================================================================
 * 3. Count 4-cycles in LDPC H matrix (structure analysis)
 * ----------------------------------------------------------------
 * 4-cycle = two variable nodes sharing ≥2 check nodes.
 * This function counts all such cycles.
 * ================================================================ */
static int factorial(int n) {
  int i, f = 1;
  for (i = 1; i <= n; i++)
    f *= i;
  return f;
}

int count_floop(int **H, int N, int wc, int wr) {
  int i, j, k, l;
  int M = (N * wc) / wr;

  int **var_nodes = (int **)malloc(N * sizeof(int *));
  for (i = 0; i < N; i++)
    var_nodes[i] = (int *)malloc(wc * sizeof(int));

  /* Construct variable-node adjacency list */
  for (j = 0; j < N; j++) {
    int idx = 0;
    for (i = 0; i < M; i++) {
      if (H[i][j]) {
        var_nodes[j][idx++] = i;
        if (idx == wc)
          break;
      }
    }
  }

  int floop = 0;

  /* Count pairs of columns that share check nodes */
  for (i = 0; i < N - 1; i++) {
    for (j = i + 1; j < N; j++) {
      int shared = 0;
      for (k = 0; k < wc; k++)
        for (l = 0; l < wc; l++)
          if (var_nodes[i][k] == var_nodes[j][l])
            shared++;

      if (shared >= 2)
        floop += factorial(shared) / (2 * factorial(shared - 2));
    }
  }

  /* free */
  for (i = 0; i < N; i++)
    free(var_nodes[i]);
  free(var_nodes);

  return floop;
}
