#include "exact.h"
#include "filters.h"
#include "grfilter.h"
#include "matrices.h"

#define LEN 32
#define FIRDEG 10
#define MSDEG 3
#define MPDEG 6
#define MAXM 3
#define MSM 8
#define MEM 4

int main(int argc, char *argv[]) {

  int n_inputs;
  double buffer_in[BATCH_SIZE][LEN];
  double buffer_out_exact[BATCH_SIZE][LEN];
  double buffer_out_mat[BATCH_SIZE][LEN];
  double buffer_out_pgf[FIRDEG][BATCH_SIZE][LEN];
  double buffer_out_mupgf[MPDEG * MAXM][BATCH_SIZE][LEN];
  double buffer_out_ms[MSDEG * MSM][BATCH_SIZE][LEN];
  double buffer_out_me[MSDEG * MEM][BATCH_SIZE][LEN];
  double buffer_temp[3][LEN];

  double diff = 0;
  double acc_error_mat = 0, acc_error_pgf[10] = {0},
         acc_error_ms[MSDEG * MSM] = {0}, acc_error_me[MSDEG * MEM] = {0},
         acc_error_mupgf[MPDEG * MAXM] = {0};

  const double *pgf_coeffs_ptr[10] = {exp32_pgf1_coeffs, exp32_pgf2_coeffs,
                                      exp32_pgf3_coeffs, exp32_pgf4_coeffs,
                                      exp32_pgf5_coeffs, exp32_pgf6_coeffs,
                                      exp32_pgf7_coeffs, exp32_pgf8_coeffs,
                                      exp32_pgf9_coeffs, exp32_pgf10_coeffs};
  const double *mupgf_coeffs_ptr[(MPDEG - 1) * MAXM] = {
      exp32_mupgf_m2l1_coeffs, exp32_mupgf_m2l2_coeffs,
      exp32_mupgf_m2l3_coeffs, exp32_mupgf_m2l4_coeffs,
      exp32_mupgf_m2l5_coeffs, exp32_mupgf_m2l6_coeffs,
      exp32_mupgf_m3l1_coeffs, exp32_mupgf_m3l2_coeffs,
      exp32_mupgf_m3l3_coeffs, exp32_mupgf_m3l4_coeffs,
      exp32_mupgf_m3l5_coeffs, exp32_mupgf_m3l6_coeffs};
  const double *me_coeffs_ptr[MSDEG * MEM] = {
      exp32_mel1m1_coeffs,
      exp32_mel1m2_coeffs,
      exp32_mel1m3_coeffs,
      NULL,
      exp32_mel2m1_coeffs,
      exp32_mel2m1_coeffs,
      NULL,
      NULL,
      exp32_mel3m1_coeffs,
      NULL,
      NULL,
      NULL,
  };
  const int *me_powers_ptr[MSDEG * MEM] = {
      exp32_mel1m1_powers,
      exp32_mel1m2_powers,
      exp32_mel1m3_powers,
      NULL,
      exp32_mel2m1_powers,
      exp32_mel2m1_powers,
      NULL,
      NULL,
      exp32_mel3m1_powers,
      NULL,
      NULL,
      NULL,
  };
  const double *ms_coeffs_ptr[24] = {
      exp32_msl1m1_coeffs, exp32_msl1m2_coeffs, exp32_msl1m3_coeffs,
      exp32_msl1m4_coeffs, exp32_msl1m5_coeffs, exp32_msl1m6_coeffs,
      exp32_msl1m7_coeffs, exp32_msl1m8_coeffs, exp32_msl2m1_coeffs,
      exp32_msl2m2_coeffs, exp32_msl2m3_coeffs, exp32_msl2m4_coeffs,
      exp32_msl2m5_coeffs, exp32_msl2m6_coeffs, exp32_msl2m7_coeffs,
      exp32_msl2m8_coeffs, exp32_msl3m1_coeffs, exp32_msl3m2_coeffs,
      exp32_msl3m3_coeffs, exp32_msl3m4_coeffs, exp32_msl3m5_coeffs,
      exp32_msl3m6_coeffs, exp32_msl3m7_coeffs, exp32_msl3m8_coeffs};
  const int *ms_powers_ptr[24] = {
      exp32_msl1m1_powers, exp32_msl1m2_powers, exp32_msl1m3_powers,
      exp32_msl1m4_powers, exp32_msl1m5_powers, exp32_msl1m6_powers,
      exp32_msl1m7_powers, exp32_msl1m8_powers, exp32_msl2m1_powers,
      exp32_msl2m2_powers, exp32_msl2m3_powers, exp32_msl2m4_powers,
      exp32_msl2m5_powers, exp32_msl2m6_powers, exp32_msl2m7_powers,
      exp32_msl2m8_powers, exp32_msl3m1_powers, exp32_msl3m2_powers,
      exp32_msl3m3_powers, exp32_msl3m4_powers, exp32_msl3m5_powers,
      exp32_msl3m6_powers, exp32_msl3m7_powers, exp32_msl3m8_powers};

  // read inputs
  FILE *fp_in = fopen(argv[1], "r");
  FILE *fp_out = fopen(argv[2], "w+");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  int n_batches = ceil((double)n_inputs / (double)BATCH_SIZE);
  int cur_batch_size = 0;
  clock_t t_temp = 0, t_exact = 0, t_mat = 0, t_pgf[10] = {0},
          t_ms[MSDEG * MSM] = {0}, t_me[MSDEG * MEM] = {0}, t_mupgf[MPDEG * 2];

  int ind = 0;

  for (int b = 0; b < n_batches; b++) {
    cur_batch_size = b < n_batches - 1 ? BATCH_SIZE : n_inputs - b * BATCH_SIZE;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < LEN; j++)
        fscanf(fp_in, "%lf", &buffer_in[i][j]);
    }

    memset(buffer_out_exact, 0, BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_mat, 0, BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_pgf, 0, FIRDEG * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_ms, 0, MSDEG * MSM * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_me, 0, MSDEG * MEM * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_out_mupgf, 0, MPDEG * 2 * BATCH_SIZE * LEN * sizeof(double));
    memset(buffer_temp, 0, MAXM * LEN * sizeof(double));

    // Exact filter
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      exact_filter_32(buffer_in[i], buffer_out_exact[i], h32_exp);
    t_exact += clock() - t_temp;

    // matrix multiplication
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      mat_times_vec(buffer_in[i], buffer_out_mat[i], exp32, LEN);
    t_mat += clock() - t_temp;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < LEN; j++) {
        diff = buffer_out_exact[i][j] - buffer_out_mat[i][j];
        acc_error_mat += diff * diff;
      }
    }

    // Polynomial graph filter (PGF)
    for (int ord = 1; ord <= FIRDEG; ord++) {
      t_temp = clock();
      for (int i = 0; i < cur_batch_size; i++)
        pgf(buffer_in[i], buffer_out_pgf[ord - 1][i], LEN, ord,
            pgf_coeffs_ptr[ord - 1], NE_LD32, MEV_LD32, Ld32_a, Ld32_w);
      t_pgf[ord - 1] += clock() - t_temp;
      for (int i = 0; i < cur_batch_size; i++) {
        for (int j = 0; j < LEN; j++) {
          diff = buffer_out_exact[i][j] - buffer_out_pgf[ord - 1][i][j];
          acc_error_pgf[ord - 1] += diff * diff;
        }
      }
    }

    // Multi-polynomial graph filter
    for (int m = 2; m <= MAXM; m++) {
      for (int l = 1; l <= MPDEG; l++) {
        memset(buffer_temp, 0, MAXM * LEN * sizeof(double));
        ind = (m - 2) * MPDEG + l - 1;
        t_temp = clock();
        for (int i = 0; i < cur_batch_size; i++) {
          pgf(buffer_in[i], buffer_temp[0], LEN, l, &mupgf_coeffs_ptr[ind][0],
              NE_LD32, 0, Ld32_a, Ld32_w);
          if (m >= 2)
            pgf_s(buffer_in[i], buffer_temp[1], LEN, l,
                  &mupgf_coeffs_ptr[ind][l + 1], NE_BD32_3, 0, Bd32_3_a,
                  Bd32_3_w);
          if (m >= 3)
            pgf_s(buffer_in[i], buffer_temp[2], LEN, l,
                  &mupgf_coeffs_ptr[ind][2 * l + 2], NE_BD32_4, 0, Bd32_4_a,
                  Bd32_4_w);
          buffer_add(buffer_temp, buffer_out_mupgf[ind][i], m, LEN);
        }
        t_mupgf[ind] += clock() - t_temp;
        for (int i = 0; i < cur_batch_size; i++) {
          for (int j = 0; j < LEN; j++) {
            diff = buffer_out_exact[i][j] - buffer_out_mupgf[ind][i][j];
            acc_error_mupgf[ind] += diff * diff;
          }
        }
      }
    }

    // MPGF--exhaustive search
    for (int l = 1; l <= MSDEG; l++) {
      for (int m = 1; m <= MEM; m++) {
        ind = (l - 1) * MEM + m - 1;
        if (!me_coeffs_ptr[ind])
          continue;
        // parse the power list
        int idx_list[MSDEG * MEM] = {0}, pow_list[MSDEG * MEM] = {0};
        get_mpgf_terms(me_powers_ptr[ind], l, m, NOPS_LD32, idx_list, pow_list);
        t_temp = clock();
        for (int i = 0; i < cur_batch_size; i++)
          mpgf(buffer_in[i], buffer_out_me[ind][i], LEN, l, m,
               me_coeffs_ptr[ind], idx_list, pow_list, nes_bd32, alists_bd32,
               wlists_bd32);
        t_me[ind] += clock() - t_temp;
        for (int i = 0; i < cur_batch_size; i++) {
          for (int j = 0; j < LEN; j++) {
            diff = buffer_out_exact[i][j] - buffer_out_me[ind][i][j];
            acc_error_me[ind] += diff * diff;
          }
        }
      }
    }

    // MPGF--OMP
    for (int l = 1; l <= MSDEG; l++) {
      for (int m = 1; m <= MSM; m++) {
        // parse the power list
        int idx_list[MSDEG * MSM] = {0}, pow_list[MSDEG * MSM] = {0};
        ind = (l - 1) * MSM + m - 1;
        get_mpgf_terms(ms_powers_ptr[ind], l, m, NOPS_LD32, idx_list, pow_list);
        t_temp = clock();
        for (int i = 0; i < cur_batch_size; i++)
          mpgf(buffer_in[i], buffer_out_ms[ind][i], LEN, l, m,
               ms_coeffs_ptr[ind], idx_list, pow_list, nes_bd32, alists_bd32,
               wlists_bd32);
        t_ms[ind] += clock() - t_temp;
        for (int i = 0; i < cur_batch_size; i++) {
          for (int j = 0; j < LEN; j++) {
            diff = buffer_out_exact[i][j] - buffer_out_ms[ind][i][j];
            acc_error_ms[ind] += diff * diff;
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
      fprintf(fp_out, "FIR order 1: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_pgf[0][i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "FIR order 2: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_pgf[1][i][j]);
      fprintf(fp_out, "\n");
      fprintf(fp_out, "Multishifts order 3 nterms 8: ");
      for (int j = 0; j < LEN; j++)
        fprintf(fp_out, "%.8lf ", buffer_out_MSM3[0][i][j]);
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
  for (int ord = 1; ord <= FIRDEG; ord++) {
    double time_pgf = ((double)t_pgf[ord - 1]) / CLOCKS_PER_SEC;
    fprintf(fp_out, "FIR filter (order = %d):    %.8lf", ord, time_pgf);
    fprintf(fp_out, " (error = %.8lf)\n",
            acc_error_pgf[ord - 1] / ((double)n_inputs));
  }
  for (int m = 2; m <= MAXM; m++) {
    for (int l = 1; l <= MPDEG; l++) {
      double time_mupgf =
          ((double)t_mupgf[(m - 2) * MPDEG + l - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "Multi-polynomial, (m = %d, order = %d):    %.8lf", m, l,
              time_mupgf);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_mupgf[(m - 2) * MSDEG + l - 1] / ((double)n_inputs));
    }
  }
  for (int l = 1; l <= MSDEG; l++) {
    for (int m = 1; m <= MEM; m++) {
      if (!me_coeffs_ptr[(l - 1) * MEM + m - 1])
        continue;
      double time_me = ((double)t_me[(l - 1) * MEM + m - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out,
              "Multishift filter, exhaustive (order = %d, m = %d):    %.8lf", l,
              m, time_me);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_me[(l - 1) * MSM + m - 1] / ((double)n_inputs));
    }
  }
  for (int l = 1; l <= MSDEG; l++) {
    for (int m = 1; m <= MSM; m++) {
      double time_ms = ((double)t_ms[(l - 1) * MSM + m - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "Multishift filter, OMP (order = %d, m = %d):    %.8lf",
              l, m, time_ms);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_ms[(l - 1) * MSM + m - 1] / ((double)n_inputs));
    }
  }

  fclose(fp_in);
  fclose(fp_out);
  return 0;
}
