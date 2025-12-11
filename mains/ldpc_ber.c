/**
 * @file ldpc_ber.c
 * @brief LDPC BER simulation over AWGN (BPSK) with SPA decoding.
 *
 * This program evaluates the information-bit BER performance of a
 * systematic LDPC code under BPSK modulation over AWGN, using the
 * Sum-Product Algorithm (SPA) decoder.
 */

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "ldpc_decoder.h"
#include "ldpc_encoder.h"

#define PI 3.141592653589793

/* ============================================================
 * Simulation parameters
 * ============================================================ */
const int N_trials = 10; /* Monte Carlo trials */
const double EbN0_min = -2.0;
const double EbN0_max = 10.0;
const double EbN0_step = 0.5;
const int max_iter_spa = 40; /* SPA maximum iteration */

/* ============================================================
 * Gaussian noise generator
 * ============================================================ */
static double rand_uniform(void) { return (rand() + 1.0) / (RAND_MAX + 2.0); }
static double randn(void) {
  double u1 = rand_uniform();
  double u2 = rand_uniform();
  return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ============================================================
 * BPSK theoretical BER
 * ============================================================ */
static double bpsk_ber(double EbN0_linear) {
  return 0.5 * erfc(sqrt(EbN0_linear));
}

/* ============================================================
 * Load matrix from CSV
 * ============================================================ */
static int load_matrix(int **mat, int rows, int cols, const char *path) {
  FILE *fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: cannot open %s\n", path);
    return -1;
  }

  char line[20000];
  for (int r = 0; r < rows; r++) {
    if (!fgets(line, sizeof(line), fp)) {
      fprintf(stderr, "ERROR: insufficient rows in %s\n", path);
      fclose(fp);
      return -2;
    }
    for (int c = 0; c < cols; c++)
      mat[r][c] = (line[c] == '1') ? 1 : 0;
  }
  fclose(fp);
  return 0;
}

/* ============================================================
 * Matrix alloc/free
 * ============================================================ */
static int **alloc_matrix_int(int rows, int cols) {
  int **m = (int **)malloc(rows * sizeof(int *));
  for (int i = 0; i < rows; i++)
    m[i] = (int *)malloc(cols * sizeof(int));
  return m;
}

static void free_matrix_int(int **m, int rows) {
  for (int i = 0; i < rows; i++)
    free(m[i]);
  free(m);
}

/* ============================================================
 * List LDPC folders under matrices/
 * ============================================================ */
static void select_ldpc_folder(char *selected_path, size_t bufsize) {
  const char *root = "matrices";
  DIR *dp = opendir(root);
  if (!dp) {
    fprintf(stderr, "ERROR: '%s' directory not found.\n", root);
    exit(1);
  }

  char *folder_list[256];
  int folder_count = 0;
  struct dirent *entry;

  printf("Searching LDPC matrix folders under '%s'...\n\n", root);

  while ((entry = readdir(dp)) != NULL) {
    if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, ".."))
      continue;

    folder_list[folder_count] = strdup(entry->d_name);
    folder_count++;
  }
  closedir(dp);

  printf("Available LDPC matrix folders:\n\n");
  for (int i = 0; i < folder_count; i++)
    printf(" [%d] %s/%s\n", i, root, folder_list[i]);

  printf("\nSelect folder index: ");
  int idx = 0;
  scanf("%d", &idx);

  snprintf(selected_path, bufsize, "%s/%s", root, folder_list[idx]);

  for (int i = 0; i < folder_count; i++)
    free(folder_list[i]);

  printf("\nUsing LDPC folder: %s\n\n", selected_path);
}

/* ============================================================
 * MAIN
 * ============================================================ */
int main(void) {
  printf("==============================================\n");
  printf("          LDPC BER Simulation (AWGN)          \n");
  printf("==============================================\n\n");

  /* 1. Select folder */
  char folder[256];
  select_ldpc_folder(folder, sizeof(folder));

  /* 2. Parse parameters from folder name */
  int N = 0, wc = 0, wr = 0;
  if (sscanf(folder, "matrices/N%d_wc%d_wr%d", &N, &wc, &wr) != 3) {
    fprintf(stderr,
            "Folder name format error. Expected matrices/N{N}_wc{wc}_wr{wr}\n");
    return 1;
  }

  int M = (N * wc) / wr;
  int K = N - M;

  printf("LDPC parameters:\n");
  printf("  N  = %d\n", N);
  printf("  K  = %d\n", K);
  printf("  M  = %d\n", M);
  printf("  wc = %d, wr = %d\n\n", wc, wr);

  /* 3. Load H,G matrices */
  int **H = alloc_matrix_int(M, N);
  int **G = alloc_matrix_int(K, N);

  char path_H[512], path_G[512];
  snprintf(path_H, sizeof(path_H), "%s/H.csv", folder);
  snprintf(path_G, sizeof(path_G), "%s/G.csv", folder);

  if (load_matrix(H, M, N, path_H) || load_matrix(G, K, N, path_G)) {
    fprintf(stderr, "Matrix load failed.\n");
    return 1;
  }

  /* 4. Create results directory */
#ifdef _WIN32
  _mkdir("results");
#else
  mkdir("results", 0777);
#endif

  /* =============================================
   * NEW: include N, wc, wr, max_iter_spa in file name
   * ============================================= */
  char csv_path[256];
  snprintf(csv_path, sizeof(csv_path),
           "results/ldpc_ber_N%d_wc%d_wr%d_iter%d_data.csv", N, wc, wr,
           max_iter_spa);

  FILE *fp = fopen(csv_path, "w");
  if (!fp) {
    fprintf(stderr, "Cannot open %s\n", csv_path);
    return 1;
  }

  fprintf(fp, "EbN0_dB,BER_info,BER_bpsk\n");

  printf("Saving results to: %s\n\n", csv_path);

  srand((unsigned int)time(NULL));

  /* buffers */
  int *inf = malloc(K * sizeof(int));
  int *code = malloc(N * sizeof(int));
  double *tx = malloc(N * sizeof(double));
  double *rx = malloc(N * sizeof(double));
  double *LLR = malloc(N * sizeof(double));
  int *ecc_hat = malloc(N * sizeof(int));
  int *inf_hat = malloc(K * sizeof(int));

  printf("EbN0_dB, BER_info, BER_bpsk\n");

  /* 6. SNR loop */
  for (double EbN0_dB = EbN0_min; EbN0_dB <= EbN0_max; EbN0_dB += EbN0_step) {

    double EbN0 = pow(10.0, EbN0_dB / 10.0);
    double R = (double)K / N;
    double sigma2 = 1.0 / (2.0 * R * EbN0);

    long long err_info = 0;
    long long total_info_bits = (long long)N_trials * K;

    for (int f = 0; f < N_trials; f++) {

      for (int i = 0; i < K; i++)
        inf[i] = rand() & 1;

      ldpc_encode(code, inf, G, N, K);

      for (int i = 0; i < N; i++)
        tx[i] = (code[i] == 1) ? +1.0 : -1.0;

      for (int i = 0; i < N; i++)
        rx[i] = tx[i] + sqrt(sigma2) * randn();

      for (int i = 0; i < N; i++)
        LLR[i] = 2.0 * rx[i] / sigma2;

      ldpc_decode_spa(LLR, ecc_hat, inf_hat, H, M, N, K, max_iter_spa);

      for (int i = 0; i < K; i++)
        if (inf[i] != inf_hat[i])
          err_info++;
    }

    double BER_info = (double)err_info / total_info_bits;
    double BER_bpsk = bpsk_ber(EbN0);

    printf("%.1f, %.10e, %.10e\n", EbN0_dB, BER_info, BER_bpsk);
    fprintf(fp, "%.1f,%.10e,%.10e\n", EbN0_dB, BER_info, BER_bpsk);
  }

  fclose(fp);

  free(inf);
  free(code);
  free(tx);
  free(rx);
  free(LLR);
  free(ecc_hat);
  free(inf_hat);

  free_matrix_int(H, M);
  free_matrix_int(G, K);

  printf("\nResults saved to %s\n", csv_path);
  return 0;
}
