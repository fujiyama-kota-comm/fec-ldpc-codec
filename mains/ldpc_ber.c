/**
 * @file ldpc_ber.c
 * @brief LDPC BER simulation over AWGN (BPSK) with SPA decoding.
 *
 * Output:
 *   results/ldpc_ber_data.csv
 *   Columns: EbN0_dB, BER_info, BER_bpsk
 *
 * LDPC matrices are selected interactively from:
 *   matrices/<folder>/H.csv
 *   matrices/<folder>/G.csv
 * where <folder> is typically N{N}_wc{wc}_wr{wr}.
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

const int N_trials = 1000;
const double EbN0_min = -2.0;
const double EbN0_max = 10.0;
const double EbN0_step = 0.5;
const int max_iter_spa = 40;

/* ============================================================
 * Utility: Gaussian random (Box–Muller)
 * ============================================================ */
static double rand_uniform(void) { return (rand() + 1.0) / (RAND_MAX + 2.0); }

static double randn(void) {
  double u1 = rand_uniform();
  double u2 = rand_uniform();
  return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ============================================================
 * BPSK theoretical BER
 *   BER = 0.5 * erfc( sqrt(Eb/N0) )
 * ============================================================ */
static double bpsk_ber(double EbN0_linear) {
  return 0.5 * erfc(sqrt(EbN0_linear));
}

/* ============================================================
 * Load 0/1 matrix from CSV (each line: "010101...")
 * ============================================================ */
static int load_matrix(int **mat, int rows, int cols, const char *path) {
  FILE *fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: cannot open %s\n", path);
    return -1;
  }
  char line[16384];
  for (int r = 0; r < rows; r++) {
    if (!fgets(line, sizeof(line), fp)) {
      fprintf(stderr, "ERROR: insufficient rows in %s\n", path);
      fclose(fp);
      return -2;
    }
    for (int c = 0; c < cols; c++) {
      mat[r][c] = (line[c] == '1') ? 1 : 0;
    }
  }
  fclose(fp);
  return 0;
}

/* ============================================================
 * Utility: allocate/free matrix
 * ============================================================ */
static int **alloc_matrix_int(int rows, int cols) {
  int **m = (int **)malloc(rows * sizeof(int *));
  if (!m) {
    fprintf(stderr, "alloc_matrix_int: malloc failed (rows)\n");
    exit(1);
  }
  for (int i = 0; i < rows; i++) {
    m[i] = (int *)malloc(cols * sizeof(int));
    if (!m[i]) {
      fprintf(stderr, "alloc_matrix_int: malloc failed (row %d)\n", i);
      exit(1);
    }
  }
  return m;
}

static void free_matrix_int(int **m, int rows) {
  for (int i = 0; i < rows; i++)
    free(m[i]);
  free(m);
}

/* ------------------------------------------------------------
 * List subdirectories under matrices/ and let user select one
 *   (i.e., matrices/<name>)
 * ------------------------------------------------------------ */
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
    /* skip "." and ".." */
    if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
      continue;

    /* NOTE: d_type may be DT_UNKNOWN on some systems,
       but for simplicity we treat all non-dot entries as candidates. */
    folder_list[folder_count] = strdup(entry->d_name);
    if (!folder_list[folder_count]) {
      fprintf(stderr, "strdup failed\n");
      exit(1);
    }
    folder_count++;
    if (folder_count >= (int)(sizeof(folder_list) / sizeof(folder_list[0])))
      break;
  }
  closedir(dp);

  if (folder_count == 0) {
    fprintf(stderr, "No LDPC folders found under '%s'.\n", root);
    exit(1);
  }

  printf("Available LDPC matrix folders:\n\n");
  for (int i = 0; i < folder_count; i++) {
    printf("  [%d] %s/%s\n", i, root, folder_list[i]);
  }
  printf("\nSelect folder index: ");

  int idx = 0;
  if (scanf("%d", &idx) != 1) {
    fprintf(stderr, "Invalid input.\n");
    exit(1);
  }

  if (idx < 0 || idx >= folder_count) {
    fprintf(stderr, "Invalid selection.\n");
    exit(1);
  }

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

  /* ------------------------------------------------------------
   * 1. Select LDPC matrix folder from ./matrices/
   * ------------------------------------------------------------ */
  char folder[256];
  select_ldpc_folder(folder, sizeof(folder));

  /* ------------------------------------------------------------
   * 2. Parse N, wc, wr from folder name (e.g., matrices/N1024_wc3_wr6)
   * ------------------------------------------------------------ */
  int N = 0, wc = 0, wr = 0;
  if (sscanf(folder, "matrices/N%d_wc%d_wr%d", &N, &wc, &wr) != 3) {
    fprintf(stderr,
            "Folder name format error. Expected matrices/N{N}_wc{wc}_wr{wr}\n");
    return 1;
  }

  int M = (N * wc) / wr;
  int K = N - M;

  printf("LDPC parameters:\n");
  printf("  N = %d\n", N);
  printf("  K = %d\n", K);
  printf("  M = %d\n", M);
  printf("  wc = %d, wr = %d\n\n", wc, wr);

  /* ------------------------------------------------------------
   * 3. Allocate matrices & load H,G
   * ------------------------------------------------------------ */
  int **H = alloc_matrix_int(M, N);
  int **G = alloc_matrix_int(K, N);

  char path_H[256], path_G[256];
  snprintf(path_H, sizeof(path_H), "%s/H.csv", folder);
  snprintf(path_G, sizeof(path_G), "%s/G.csv", folder);

  printf("Loading: %s\n", path_H);
  printf("Loading: %s\n", path_G);

  if (load_matrix(H, M, N, path_H) || load_matrix(G, K, N, path_G)) {
    fprintf(stderr, "Matrix load failed.\n");
    return 1;
  }

  printf("Matrix load OK.\n\n");

  /* ------------------------------------------------------------
   * 4. Create results directory
   * ------------------------------------------------------------ */
#ifdef _WIN32
  _mkdir("results");
#else
  mkdir("results", 0777);
#endif

  FILE *fp = fopen("results/ldpc_ber_data.csv", "w");
  if (!fp) {
    fprintf(stderr, "Cannot open results/ldpc_ber_data.csv\n");
    return 1;
  }
  /* Only Info BER & uncoded BPSK */
  fprintf(fp, "EbN0_dB,BER_info,BER_bpsk\n");

  srand((unsigned int)time(NULL));

  /* buffers */
  int *inf = (int *)malloc(K * sizeof(int));
  int *code = (int *)malloc(N * sizeof(int));
  double *tx = (double *)malloc(N * sizeof(double));
  double *rx = (double *)malloc(N * sizeof(double));
  double *LLR = (double *)malloc(N * sizeof(double));
  int *ecc_hat = (int *)malloc(N * sizeof(int));
  int *inf_hat = (int *)malloc(K * sizeof(int));

  if (!inf || !code || !tx || !rx || !LLR || !ecc_hat || !inf_hat) {
    fprintf(stderr, "Memory allocation failed.\n");
    return 1;
  }

  printf("EbN0_dB, BER_info, BER_bpsk\n");

  /* ------------------------------------------------------------
   * 6. SNR loop
   * ------------------------------------------------------------ */
  for (double EbN0_dB = EbN0_min; EbN0_dB <= EbN0_max; EbN0_dB += EbN0_step) {

    double EbN0 = pow(10.0, EbN0_dB / 10.0);
    double R = (double)K / N;
    double sigma2 = 1.0 / (2.0 * R * EbN0);
    double sigma = sqrt(sigma2);

    long long err_info = 0;
    long long total_info_bits = (long long)N_trials * K;

    for (int f = 0; f < N_trials; f++) {

      /* info bits */
      for (int i = 0; i < K; i++)
        inf[i] = rand() & 1;

      /* encode */
      ldpc_encode(code, inf, G, N, K);

      /* BPSK mapping */
      for (int i = 0; i < N; i++)
        tx[i] = (code[i] == 1) ? +1.0 : -1.0;

      /* AWGN */
      for (int i = 0; i < N; i++)
        rx[i] = tx[i] + sigma * randn();

      /* LLR for BPSK over AWGN: LLR = 2y / sigma^2 */
      for (int i = 0; i < N; i++)
        LLR[i] = 2.0 * rx[i] / sigma2;

      /* SPA decoding */
      ldpc_decode_spa(LLR, ecc_hat, inf_hat, H, M, N, K, max_iter_spa);

      /* count info bit errors */
      for (int i = 0; i < K; i++)
        if (inf[i] != inf_hat[i])
          err_info++;
    }

    double BER_info = (double)err_info / (double)total_info_bits;
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

  printf("\nResults saved to results/ldpc_ber_data.csv\n");
  return 0;
}
