/**
 * @file gene_hg.c
 * @brief LDPC H/G matrix generator (Gallager construction + G from H)
 *
 * This tool:
 *   1. Generates LDPC parity-check matrix H via Gallager's method
 *   2. Constructs systematic generator matrix G from H
 *   3. Evaluates 4-cycles in H
 *   4. Searches for the smallest-4-cycle H/G pair
 *   5. Saves the best matrices into CSV files
 *
 * NOTE:
 *   This code repeatedly generates LDPC matrices and tracks the structure with
 *   the minimum number of 4-cycles. The algorithm is computationally heavy for
 *   large N. Adjust loop_count_max as needed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> /* mkdir() for Linux / Windows (mingw) */
#include <sys/types.h>
#include <time.h>


#include "ldpc_matrix.h"

/* portable mkdir wrapper */
static void make_dir(const char *d) {
#ifdef _WIN32
  _mkdir(d);
#else
  mkdir(d, 0755);
#endif
}

int main(void) {
  srand((unsigned int)time(NULL));

  printf("==============================================\n");
  printf("       LDPC Matrix Generator (Gallager)       \n");
  printf("==============================================\n\n");

  /* ----------------------------------------------------------------------
   * User input
   * ---------------------------------------------------------------------- */
  int N, wc, wr;

  printf("Codeword length N: ");
  scanf("%d", &N);

  printf("Column weight wc (small: 2 or 3): ");
  scanf("%d", &wc);

  printf("Row weight wr (larger than wc): ");
  scanf("%d", &wr);

  int M = (N * wc) / wr; /* Number of parity equations */
  int K = N - M;         /* Number of information bits */
  double R = (double)K / N;

  printf("\nRate R = %.5f (K = %d, M = %d)\n\n", R, K, M);

  /* ----------------------------------------------------------------------
   * Make output directory
   * ---------------------------------------------------------------------- */
  char dirpath[128];
  sprintf(dirpath, "matrices/N%d_wc%d_wr%d", N, wc, wr);

  make_dir("matrices");
  make_dir(dirpath);

  /* Output file paths */
  char path_H[256], path_G[256], path_info[256];
  sprintf(path_H, "%s/H.csv", dirpath);
  sprintf(path_G, "%s/G.csv", dirpath);
  sprintf(path_info, "%s/info.txt", dirpath);

  /* ----------------------------------------------------------------------
   * Allocate matrices
   * ---------------------------------------------------------------------- */
  int **H = (int **)malloc(M * sizeof(int *));
  int **H_best = (int **)malloc(M * sizeof(int *));
  for (int i = 0; i < M; i++) {
    H[i] = (int *)malloc(N * sizeof(int));
    H_best[i] = (int *)malloc(N * sizeof(int));
  }

  int **G = (int **)malloc(K * sizeof(int *));
  int **G_best = (int **)malloc(K * sizeof(int *));
  for (int i = 0; i < K; i++) {
    G[i] = (int *)malloc(N * sizeof(int));
    G_best[i] = (int *)malloc(N * sizeof(int));
  }

  /* ----------------------------------------------------------------------
   * Search best H/G matrices (min 4-cycles)
   * ---------------------------------------------------------------------- */
  const int loop_count_max = 2000000000; /* WARNING: extremely heavy */
  const double print_interval_sec = 1.0;

  int loop;
  int best_floop = -1;
  long long floop_sum = 0;

  clock_t t_start = clock();
  clock_t t_last_print = clock();

  printf("Searching for best H/G matrices (min 4-cycles)...\n");

  for (loop = 1; loop <= loop_count_max; loop++) {

    /* Generate H and G */
    generate_Hmatrix(H, N, wc, wr);
    generate_Gmatrix(H, G, N, wc, wr);

    /* Count 4-cycles */
    int floop = count_floop(H, N, wc, wr);
    floop_sum += floop;

    /* update best result */
    if (best_floop == -1 || floop < best_floop) {
      best_floop = floop;

      /* copy best matrices */
      for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
          H_best[i][j] = H[i][j];

      for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
          G_best[i][j] = G[i][j];
    }

    /* periodic saving */
    clock_t t_now = clock();
    double elapsed = (double)(t_now - t_last_print) / CLOCKS_PER_SEC;

    if (loop == 1 || elapsed > print_interval_sec) {

      t_last_print = t_now;

      /* write H */
      FILE *fp = fopen(path_H, "w");
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++)
          fprintf(fp, "%d", H_best[i][j]);
        fprintf(fp, "\n");
      }
      fclose(fp);

      /* write G */
      fp = fopen(path_G, "w");
      for (int i = 0; i < K; i++) {
        for (int j = 0; j < N; j++)
          fprintf(fp, "%d", G_best[i][j]);
        fprintf(fp, "\n");
      }
      fclose(fp);

      /* write info */
      fp = fopen(path_info, "w");
      fprintf(fp, "LDPC Matrix Generation Status\n");
      fprintf(fp, "Code rate R = %.5f\n", R);
      fprintf(fp, "N = %d\n", N);
      fprintf(fp, "wc = %d\n", wc);
      fprintf(fp, "wr = %d\n", wr);
      fprintf(fp, "Loop count = %d\n", loop);
      fprintf(fp, "Best 4-cycles = %d\n", best_floop);
      fprintf(fp, "Average 4-cycles = %.3f\n", (double)floop_sum / loop);
      fclose(fp);

      printf("[Loop %d] Best 4-cycles = %d, Avg = %.3f\n", loop, best_floop,
             (double)floop_sum / loop);
    }
  }

  /* ----------------------------------------------------------------------
   * Free memory
   * ---------------------------------------------------------------------- */
  for (int i = 0; i < M; i++) {
    free(H[i]);
    free(H_best[i]);
  }
  free(H);
  free(H_best);

  for (int i = 0; i < K; i++) {
    free(G[i]);
    free(G_best[i]);
  }
  free(G);
  free(G_best);

  printf("\nGeneration completed.\n");
  printf("Files saved under directory: %s\n", dirpath);

  return 0;
}
