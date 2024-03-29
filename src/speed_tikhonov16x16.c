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

  double diff;
  double acc_error_mat = 0, acc_error_pgf[PGFDEG] = {0},
         acc_error_cpgf[PGFDEG] = {0}, acc_error_mpgf[MPGFDEG * MAXM] = {0},
         acc_error_arma[ABDEGS * NITS] = {0};

  const double *pgf_coeffs_ptr[PGFDEG] = {
      tik16x16_pgf1_coeffs, tik16x16_pgf2_coeffs, tik16x16_pgf3_coeffs,
      tik16x16_pgf4_coeffs, tik16x16_pgf5_coeffs, tik16x16_pgf6_coeffs,
      tik16x16_pgf7_coeffs, tik16x16_pgf8_coeffs, tik16x16_pgf9_coeffs,
      tik16x16_pgf10_coeffs};
  const double *cpgf_coeffs_ptr[PGFDEG] = {
      tik16x16_cpgf1_coeffs, tik16x16_cpgf2_coeffs, tik16x16_cpgf3_coeffs,
      tik16x16_cpgf4_coeffs, tik16x16_cpgf5_coeffs, tik16x16_cpgf6_coeffs,
      tik16x16_cpgf7_coeffs, tik16x16_cpgf8_coeffs, tik16x16_cpgf9_coeffs,
      tik16x16_cpgf10_coeffs};
  const double *mpgf_coeffs_ptr[MPGFDEG * MAXM] = {
      tik16x16_mpgfl1m1_coeffs, tik16x16_mpgfl1m2_coeffs, tik16x16_mpgfl1m3_coeffs,
      tik16x16_mpgfl1m4_coeffs, tik16x16_mpgfl1m5_coeffs, tik16x16_mpgfl1m6_coeffs,
      tik16x16_mpgfl1m7_coeffs, tik16x16_mpgfl1m8_coeffs, tik16x16_mpgfl2m1_coeffs,
      tik16x16_mpgfl2m2_coeffs, tik16x16_mpgfl2m3_coeffs, tik16x16_mpgfl2m4_coeffs,
      tik16x16_mpgfl2m5_coeffs, tik16x16_mpgfl2m6_coeffs, tik16x16_mpgfl2m7_coeffs,
      tik16x16_mpgfl2m8_coeffs};
  const int *mpgf_powers_ptr[MPGFDEG * MAXM] = {
      tik16x16_mpgfl1m1_powers, tik16x16_mpgfl1m2_powers, tik16x16_mpgfl1m3_powers,
      tik16x16_mpgfl1m4_powers, tik16x16_mpgfl1m5_powers, tik16x16_mpgfl1m6_powers,
      tik16x16_mpgfl1m7_powers, tik16x16_mpgfl1m8_powers, tik16x16_mpgfl2m1_powers,
      tik16x16_mpgfl2m2_powers, tik16x16_mpgfl2m3_powers, tik16x16_mpgfl2m4_powers,
      tik16x16_mpgfl2m5_powers, tik16x16_mpgfl2m6_powers, tik16x16_mpgfl2m7_powers,
      tik16x16_mpgfl2m8_powers};
  const double *arma_b_ptr[ABDEGS] = {tik16x16_armad1_b, tik16x16_armad2_b,
                               tik16x16_armad3_b};
  const double *arma_a_ptr[ABDEGS] = {tik16x16_armad1_a, tik16x16_armad2_a,
                               tik16x16_armad3_a};

  // read inputs
  FILE *fp_in = fopen(argv[1], "r");
  FILE *fp_out = fopen(argv[2], "w+");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  int n_batches = ceil((double)n_inputs / (double)BATCH_SIZE);
  int cur_batch_size = 0;
  clock_t t_temp = 0, t_exact = 0, t_mat = 0, t_pgf[PGFDEG] = {0},
          t_cpgf[PGFDEG] = {0}, t_mpgf[MPGFDEG * MAXM] = {0},
          t_arma[ABDEGS * NITS] = {0};

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
    memset(buffer_out_cpgf, 0, PGFDEG * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_mpgf, 0, MPGFDEG * MAXM * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_arma, 0,
           ABDEGS * NITS * BATCH_SIZE * LEN * sizeof(double));

    // Exact filter
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      exact_filter_16x16(buffer_in[i], buffer_out_exact[i], h16x16_tik);
    t_exact += clock() - t_temp;

    // matrix multiplication
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      mat_times_vec(buffer_in[i], buffer_out_mat[i], tik16x16, LEN);
    t_mat += clock() - t_temp;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < LEN; j++) {
        diff = buffer_out_exact[i][j] - buffer_out_mat[i][j];
        acc_error_mat += diff * diff;
      }
    }

    // Polynomial graph filter (PGF)
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

    // Chebyshev polynomial graph filter (CPGF)
    for (int ord = 1; ord <= PGFDEG; ord++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        chebyshev_gf(buffer_in[i], buffer_out_cpgf[ord - 1][i], LEN, ord,
                     cpgf_coeffs_ptr[ord - 1], NE_LDD16X16, MEV_LDD16X16, Ldd16x16_a,
                     Ldd16x16_w);
      t_cpgf[ord - 1] += clock() - t_temp;
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < LEN; j++) {
          diff = buffer_out_exact[i][j] - buffer_out_cpgf[ord - 1][i][j];
          acc_error_cpgf[ord - 1] += diff * diff;
        }
      }
    }

    // MPGF--OMP
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

    // ARMA, CG
    for (int d = 0; d < ABDEGS; d++) {
      for (int it = 1; it <= NITS; it++) {
        ind = d * NITS + it - 1;
        t_temp = clock();
        for (int i = 0; i < cur_batch_size; i++)
          armagf_cg(buffer_in[i], buffer_out_arma[ind][i], LEN, it, kbs[d],
                    arma_b_ptr[d], kas[d], arma_a_ptr[d], NE_LDD16X16, Ldd16x16_a,
                    Ldd16x16_w);
        t_arma[ind] += clock() - t_temp;
        for (int i = 0; i < cur_batch_size; i++) {
          for (int j = 0; j < LEN; j++) {
            diff = buffer_out_exact[i][j] - buffer_out_arma[ind][i][j];
            acc_error_arma[ind] += diff * diff;
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
    fprintf(fp_out, "PGF (order = %d):    %.8lf", ord, time_pgf);
    fprintf(fp_out, " (error = %.8lf)\n",
            acc_error_pgf[ord - 1] / ((double)n_inputs));
  }
  for (int ord = 1; ord <= PGFDEG; ord++) {
    double time_cpgf = ((double)t_cpgf[ord - 1]) / CLOCKS_PER_SEC;
    fprintf(fp_out, "CPGF (order = %d):    %.8lf", ord, time_cpgf);
    fprintf(fp_out, " (error = %.8lf)\n",
            acc_error_cpgf[ord - 1] / ((double)n_inputs));
  }
  for (int l = 1; l <= MPGFDEG; l++) {
    for (int m = 1; m <= MAXM; m++) {
      double time_mpgf = ((double)t_mpgf[(l - 1) * MAXM + m - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "MPGF, OMP (order = %d, m = %d):    %.8lf",
              l, m, time_mpgf);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_mpgf[(l - 1) * MAXM + m - 1] / ((double)n_inputs));
    }
  }
  for (int d = 0; d < ABDEGS; d++) {
    for (int it = 1; it <= NITS; it++) {
      double time_arma = ((double)t_arma[d * NITS + it - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "ARMA, CG (Q = %d, P = %d, iter = %d):    %.8lf", kbs[d],
              kas[d], it, time_arma);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_arma[d * NITS + it - 1] / ((double)n_inputs));
    }
  }

  fclose(fp_in);
  fclose(fp_out);
  return 0;
}
