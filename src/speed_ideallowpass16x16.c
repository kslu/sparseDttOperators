#include "exact.h"
#include "filters.h"
#include "grfilter.h"
#include "matrices.h"

#define LEN 256
#define PGFDEG 10
#define MPGFDEG 2
#define MAXM 8
#define ABDEGS 3
#define NITS 6

int main(int argc, char *argv[]) {

  int n_inputs;
  double buffer_in[BATCH_SIZE][LEN];
  double buffer_out_exact[BATCH_SIZE][LEN];
  double buffer_out_mat[BATCH_SIZE][LEN];
  double buffer_out_pgf[PGFDEG][BATCH_SIZE][LEN];
  double buffer_out_cpgf[PGFDEG][BATCH_SIZE][LEN];
  double buffer_out_mpgf[MPGFDEG * MAXM][BATCH_SIZE][LEN];
  double buffer_out_arma[ABDEGS * NITS][BATCH_SIZE][LEN];

  int kbs[ABDEGS] = {1, 2, 3};
  int kas[ABDEGS] = {1, 2, 3};

  double diff = 0;
  double acc_error_mat = 0, acc_error_pgf[PGFDEG] = {0},
         acc_error_pgf_mm[PGFDEG] = {0}, acc_error_mpgf[MPGFDEG * MAXM] = {0},
         acc_error_mpgf_mm[MPGFDEG * MAXM] = {0};

  const double *pgf_coeffs_ptr[PGFDEG] = {
      lp16x16_pgf1_coeffs, lp16x16_pgf2_coeffs, lp16x16_pgf3_coeffs,
      lp16x16_pgf4_coeffs, lp16x16_pgf5_coeffs, lp16x16_pgf6_coeffs,
      lp16x16_pgf7_coeffs, lp16x16_pgf8_coeffs, lp16x16_pgf9_coeffs,
      lp16x16_pgf10_coeffs};
  const double *pgf_minimax_coeffs_ptr[PGFDEG] = {
      lp_minimax16x16_pgf1_coeffs, lp_minimax16x16_pgf2_coeffs,
      lp_minimax16x16_pgf3_coeffs, lp_minimax16x16_pgf4_coeffs,
      lp_minimax16x16_pgf5_coeffs, lp_minimax16x16_pgf6_coeffs,
      lp_minimax16x16_pgf7_coeffs, lp_minimax16x16_pgf8_coeffs,
      lp_minimax16x16_pgf9_coeffs, lp_minimax16x16_pgf10_coeffs};
  const double *mpgf_coeffs_ptr[MPGFDEG * MAXM] = {
      lp16x16_mpgfl1m1_coeffs, lp16x16_mpgfl1m2_coeffs, lp16x16_mpgfl1m3_coeffs,
      lp16x16_mpgfl1m4_coeffs, lp16x16_mpgfl1m5_coeffs, lp16x16_mpgfl1m6_coeffs,
      lp16x16_mpgfl1m7_coeffs, lp16x16_mpgfl1m8_coeffs, lp16x16_mpgfl2m1_coeffs,
      lp16x16_mpgfl2m2_coeffs, lp16x16_mpgfl2m3_coeffs, lp16x16_mpgfl2m4_coeffs,
      lp16x16_mpgfl2m5_coeffs, lp16x16_mpgfl2m6_coeffs, lp16x16_mpgfl2m7_coeffs,
      lp16x16_mpgfl2m8_coeffs};
  const int *mpgf_powers_ptr[MPGFDEG * MAXM] = {
      lp16x16_mpgfl1m1_powers, lp16x16_mpgfl1m2_powers, lp16x16_mpgfl1m3_powers,
      lp16x16_mpgfl1m4_powers, lp16x16_mpgfl1m5_powers, lp16x16_mpgfl1m6_powers,
      lp16x16_mpgfl1m7_powers, lp16x16_mpgfl1m8_powers, lp16x16_mpgfl2m1_powers,
      lp16x16_mpgfl2m2_powers, lp16x16_mpgfl2m3_powers, lp16x16_mpgfl2m4_powers,
      lp16x16_mpgfl2m5_powers, lp16x16_mpgfl2m6_powers, lp16x16_mpgfl2m7_powers,
      lp16x16_mpgfl2m8_powers};
  const double *mpgf_minimax_coeffs_ptr[MPGFDEG * MAXM] = {
      lp_minimax16x16_mpgfl1m1_coeffs, lp_minimax16x16_mpgfl1m2_coeffs,
      lp_minimax16x16_mpgfl1m3_coeffs, lp_minimax16x16_mpgfl1m4_coeffs,
      lp_minimax16x16_mpgfl1m5_coeffs, lp_minimax16x16_mpgfl1m6_coeffs,
      lp_minimax16x16_mpgfl1m7_coeffs, lp_minimax16x16_mpgfl1m8_coeffs,
      lp_minimax16x16_mpgfl2m1_coeffs, lp_minimax16x16_mpgfl2m2_coeffs,
      lp_minimax16x16_mpgfl2m3_coeffs, lp_minimax16x16_mpgfl2m4_coeffs,
      lp_minimax16x16_mpgfl2m5_coeffs, lp_minimax16x16_mpgfl2m6_coeffs,
      lp_minimax16x16_mpgfl2m7_coeffs, lp_minimax16x16_mpgfl2m8_coeffs};
  const int *mpgf_minimax_powers_ptr[MPGFDEG * MAXM] = {
      lp_minimax16x16_mpgfl1m1_powers, lp_minimax16x16_mpgfl1m2_powers,
      lp_minimax16x16_mpgfl1m3_powers, lp_minimax16x16_mpgfl1m4_powers,
      lp_minimax16x16_mpgfl1m5_powers, lp_minimax16x16_mpgfl1m6_powers,
      lp_minimax16x16_mpgfl1m7_powers, lp_minimax16x16_mpgfl1m8_powers,
      lp_minimax16x16_mpgfl2m1_powers, lp_minimax16x16_mpgfl2m2_powers,
      lp_minimax16x16_mpgfl2m3_powers, lp_minimax16x16_mpgfl2m4_powers,
      lp_minimax16x16_mpgfl2m5_powers, lp_minimax16x16_mpgfl2m6_powers,
      lp_minimax16x16_mpgfl2m7_powers, lp_minimax16x16_mpgfl2m8_powers};

  // read inputs
  FILE *fp_in = fopen(argv[1], "r");
  FILE *fp_out = fopen(argv[2], "w+");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  int n_batches = ceil((double)n_inputs / (double)BATCH_SIZE);
  int cur_batch_size = 0;
  clock_t t_temp = 0, t_exact = 0, t_mat = 0, t_pgf[PGFDEG] = {0},
          t_pgf_mm[PGFDEG] = {0}, t_mpgf[MPGFDEG * MAXM] = {0},
          t_mpgf_mm[MPGFDEG * MAXM] = {0};

  int ind = 0;

  for (int b = 0; b < n_batches; b++) {
    cur_batch_size = b < n_batches - 1 ? BATCH_SIZE : n_inputs - b * BATCH_SIZE;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < LEN; j++)
        fscanf(fp_in, "%lf", &buffer_in[i][j]);
    }

    memset(buffer_out_exact, 0, BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_mat, 0, BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_pgf, 0, PGFDEG * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_mpgf, 0, MPGFDEG * MAXM * BATCH_SIZE * LEN * sizeof(double));

    // Exact filter
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      exact_filter_16x16(buffer_in[i], buffer_out_exact[i], h16x16_lp);
    t_exact += clock() - t_temp;

    // matrix multiplication
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      mat_times_vec(buffer_in[i], buffer_out_mat[i], lp16x16, LEN);
    t_mat += clock() - t_temp;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < LEN; j++) {
        diff = buffer_out_exact[i][j] - buffer_out_mat[i][j];
        acc_error_mat += diff * diff;
      }
    }

    // PGF, LS
    for (int ord = 1; ord <= PGFDEG; ord++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        pgf(buffer_in[i], buffer_out_pgf[ord - 1][i], LEN, ord,
            pgf_coeffs_ptr[ord - 1], NE_LDD16X16, MEV_LDD16X16, Ldd16x16_a, Ldd16x16_w);
      t_pgf[ord - 1] += clock() - t_temp;
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < LEN; j++) {
          diff = buffer_out_exact[i][j] - buffer_out_pgf[ord - 1][i][j];
          acc_error_pgf[ord - 1] += diff * diff;
        }
      }
    }

    memset(buffer_out_pgf, 0, PGFDEG * BATCH_SIZE * LEN * sizeof(double));

    // PGF, minimax
    for (int ord = 1; ord <= PGFDEG; ord++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        pgf(buffer_in[i], buffer_out_pgf[ord - 1][i], LEN, ord,
            pgf_minimax_coeffs_ptr[ord - 1], NE_LDD16X16, MEV_LDD16X16,
            Ldd16x16_a, Ldd16x16_w);
      t_pgf_mm[ord - 1] += clock() - t_temp;
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < LEN; j++) {
          diff = buffer_out_exact[i][j] - buffer_out_pgf[ord - 1][i][j];
          acc_error_pgf_mm[ord - 1] += diff * diff;
        }
      }
    }

    // MPGF, OMP, LS
    for (int l = 1; l <= MPGFDEG; l++) {
      for (int m = 1; m <= MAXM; m++) {
        // parse the power list
        int idx_list[MPGFDEG * MAXM] = {0}, pow_list[MPGFDEG * MAXM] = {0};
        ind = (l - 1) * MAXM + m - 1;
        get_mpgf_terms(mpgf_powers_ptr[ind], l, m, NOPS_LDD16X16, idx_list,
                       pow_list);
        t_temp = clock();
        for (int i = 0; i < cur_batch_size; i++)
          mpgf(buffer_in[i], buffer_out_mpgf[ind][i], LEN, l, m,
               mpgf_coeffs_ptr[ind], idx_list, pow_list, nes_bdd16x16,
               alists_bdd16x16, wlists_bdd16x16);
        t_mpgf[ind] += clock() - t_temp;
        for (int i = 0; i < cur_batch_size; i++) {
          for (int j = 0; j < LEN; j++) {
            diff = buffer_out_exact[i][j] - buffer_out_mpgf[ind][i][j];
            acc_error_mpgf[ind] += diff * diff;
          }
        }
      }
    }

    memset(buffer_out_mpgf, 0, MPGFDEG * MAXM * BATCH_SIZE * LEN * sizeof(double));

    // MPGF, OMP, minimax
    for (int l = 1; l <= MPGFDEG; l++) {
      for (int m = 1; m <= MAXM; m++) {
        // parse the power list
        int idx_list[MPGFDEG * MAXM] = {0}, pow_list[MPGFDEG * MAXM] = {0};
        ind = (l - 1) * MAXM + m - 1;
        get_mpgf_terms(mpgf_minimax_powers_ptr[ind], l, m, NOPS_LDD16X16,
                       idx_list, pow_list);
        t_temp = clock();
        for (int i = 0; i < cur_batch_size; i++)
          mpgf(buffer_in[i], buffer_out_mpgf[ind][i], LEN, l, m,
               mpgf_minimax_coeffs_ptr[ind], idx_list, pow_list, nes_bdd16x16,
               alists_bdd16x16, wlists_bdd16x16);
        t_mpgf_mm[ind] += clock() - t_temp;
        for (int i = 0; i < cur_batch_size; i++) {
          for (int j = 0; j < LEN; j++) {
            diff = buffer_out_exact[i][j] - buffer_out_mpgf[ind][i][j];
            acc_error_mpgf_mm[ind] += diff * diff;
          }
        }
      }
    }

#if CONFIG_DEBUG
    // write output GFT coefficients
    for (int i = 0; i < cur_batch_size; i++) {
      fprintf(fp_out, "Input #%d: \n", b * BATCH_SIZE + i);
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_in[i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "Exact filter: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_exact[i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "PGF order 1: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_pgf[0][i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "PGF order 2: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_pgf[1][i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "PGF order 3: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_pgf[2][i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "MPGF order 3 nterms 8: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_mpgf[MAXM * 2][i][j]);
      fprintf(fp_out, "\n");
    }
#endif
  }

  // write runtime
  double time_exact = ((double)t_exact) / CLOCKS_PER_SEC;
  double time_mat = ((double)t_mat) / CLOCKS_PER_SEC;
  fprintf(fp_out, "#input = %d\n", n_inputs);
  fprintf(fp_out, "Exact filter:    %.8lf\n", time_exact);
  fprintf(fp_out, "Matrix filter:    %.8lf", time_mat);
  fprintf(fp_out, " (error = %.8lf)\n", acc_error_mat / ((double)n_inputs));
  for (int ord = 1; ord <= PGFDEG; ord++) {
    double time_pgf = ((double)t_pgf[ord - 1]) / CLOCKS_PER_SEC;
    fprintf(fp_out, "PGF, LS (order = %d):    %.8lf", ord, time_pgf);
    fprintf(fp_out, " (error = %.8lf)\n",
            acc_error_pgf[ord - 1] / ((double)n_inputs));
  }
  for (int ord = 1; ord <= PGFDEG; ord++) {
    double time_pgf = ((double)t_pgf_mm[ord - 1]) / CLOCKS_PER_SEC;
    fprintf(fp_out, "PGF, minimax (order = %d):    %.8lf", ord, time_pgf);
    fprintf(fp_out, " (error = %.8lf)\n",
            acc_error_pgf_mm[ord - 1] / ((double)n_inputs));
  }
  for (int l = 1; l <= MPGFDEG; l++) {
    for (int m = 1; m <= MAXM; m++) {
      double time_mpgf = ((double)t_mpgf[(l - 1) * MAXM + m - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "MPGF, OMP, LS (order = %d, m = %d):    %.8lf",
              l, m, time_mpgf);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_mpgf[(l - 1) * MAXM + m - 1] / ((double)n_inputs));
    }
  }
  for (int l = 1; l <= MPGFDEG; l++) {
    for (int m = 1; m <= MAXM; m++) {
      double time_mpgf = ((double)t_mpgf_mm[(l - 1) * MAXM + m - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "MPGF, OMP, minimax (order = %d, m = %d):    %.8lf",
              l, m, time_mpgf);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_mpgf_mm[(l - 1) * MAXM + m - 1] / ((double)n_inputs));
    }
  }

  fclose(fp_in);
  fclose(fp_out);
  return 0;
}
