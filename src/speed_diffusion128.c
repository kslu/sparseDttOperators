#include "exact.h"
#include "filters.h"
#include "grfilter.h"
#include "matrices.h"

#define LEN 128
#define PGFDEG 10
#define MSDEG 3
#define MSM 8

int main(int argc, char *argv[]) {

  int n_inputs;
  double buffer_in[BATCH_SIZE][LEN];
  double buffer_out_exact[BATCH_SIZE][LEN];
  double buffer_out_pgf[PGFDEG][BATCH_SIZE][LEN];
  double buffer_out_ms[MSDEG * MSM][BATCH_SIZE][LEN];

  double diff = 0;
  double acc_error_pgf[10] = {0}, acc_error_ms[3 * 8] = {0};

  const double *pgf_coeffs_ptr[10] = {
      diff128_pgf1_coeffs, diff128_pgf2_coeffs, diff128_pgf3_coeffs,
      diff128_pgf4_coeffs, diff128_pgf5_coeffs, diff128_pgf6_coeffs,
      diff128_pgf7_coeffs, diff128_pgf8_coeffs, diff128_pgf9_coeffs,
      diff128_pgf10_coeffs};
  const double *ms_coeffs_ptr[24] = {
      diff128_msl1m1_coeffs, diff128_msl1m2_coeffs, diff128_msl1m3_coeffs,
      diff128_msl1m4_coeffs, diff128_msl1m5_coeffs, diff128_msl1m6_coeffs,
      diff128_msl1m7_coeffs, diff128_msl1m8_coeffs, diff128_msl2m1_coeffs,
      diff128_msl2m2_coeffs, diff128_msl2m3_coeffs, diff128_msl2m4_coeffs,
      diff128_msl2m5_coeffs, diff128_msl2m6_coeffs, diff128_msl2m7_coeffs,
      diff128_msl2m8_coeffs, diff128_msl3m1_coeffs, diff128_msl3m2_coeffs,
      diff128_msl3m3_coeffs, diff128_msl3m4_coeffs, diff128_msl3m5_coeffs,
      diff128_msl3m6_coeffs, diff128_msl3m7_coeffs, diff128_msl3m8_coeffs};
  const int *ms_powers_ptr[24] = {
      diff128_msl1m1_powers, diff128_msl1m2_powers, diff128_msl1m3_powers,
      diff128_msl1m4_powers, diff128_msl1m5_powers, diff128_msl1m6_powers,
      diff128_msl1m7_powers, diff128_msl1m8_powers, diff128_msl2m1_powers,
      diff128_msl2m2_powers, diff128_msl2m3_powers, diff128_msl2m4_powers,
      diff128_msl2m5_powers, diff128_msl2m6_powers, diff128_msl2m7_powers,
      diff128_msl2m8_powers, diff128_msl3m1_powers, diff128_msl3m2_powers,
      diff128_msl3m3_powers, diff128_msl3m4_powers, diff128_msl3m5_powers,
      diff128_msl3m6_powers, diff128_msl3m7_powers, diff128_msl3m8_powers};

  // read inputs
  FILE *fp_in = fopen(argv[1], "r");
  FILE *fp_out = fopen(argv[2], "w+");
  if (fp_in != NULL)
    fscanf(fp_in, "%d", &n_inputs);

  int n_batches = ceil((double)n_inputs / (double)BATCH_SIZE);
  int cur_batch_size = 0;
  clock_t t_temp = 0, t_mat = 0, t_pgf[10] = {0}, t_ms[3 * 8] = {0};

  int ind = 0;

  for (int b = 0; b < n_batches; b++) {
    cur_batch_size = b < n_batches - 1 ? BATCH_SIZE : n_inputs - b * BATCH_SIZE;
    for (int i = 0; i < cur_batch_size; i++) {
      for (int j = 0; j < LEN; j++)
        fscanf(fp_in, "%lf", &buffer_in[i][j]);
    }

    // matrix multiplication
    t_temp = clock();
    for (int i = 0; i < cur_batch_size; i++)
      mat_times_vec(buffer_in[i], buffer_out_exact[i], diff128, LEN);
    t_mat += clock() - t_temp;

    // Polynomial graph filter (PGF)
    for (int ord = 1; ord <= PGFDEG; ord++) {
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

    // MPGF
    for (int l = 1; l <= MSDEG; l++) {
      for (int m = 1; m <= MSM; m++) {
        // parse the power list
        ind = (l - 1) * MSM + m - 1;
        int idx_list[24] = {0}, pow_list[24] = {0};
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
  double time_mat = ((double)t_mat) / CLOCKS_PER_SEC;
  fprintf(fp_out, "#input = %d\n", n_inputs);
  fprintf(fp_out, "Exact filter:    %.8lf\n", time_mat);
  for (int ord = 1; ord <= PGFDEG; ord++) {
    double time_pgf = ((double)t_pgf[ord - 1]) / CLOCKS_PER_SEC;
    fprintf(fp_out, "FIR filter (order = %d):    %.8lf", ord, time_pgf);
    fprintf(fp_out, " (error = %.8lf)\n",
            acc_error_pgf[ord - 1] / ((double)n_inputs));
  }
  for (int l = 1; l <= MSDEG; l++) {
    for (int m = 1; m <= MSM; m++) {
      double time_ms = ((double)t_ms[(l - 1) * MSM + m - 1]) / CLOCKS_PER_SEC;
      fprintf(fp_out, "Multishift filter (order = %d, m = %d):    %.8lf", l, m,
              time_ms);
      fprintf(fp_out, " (error = %.8lf)\n",
              acc_error_ms[(l - 1) * MSM + m - 1] / ((double)n_inputs));
    }
  }

  fclose(fp_in);
  fclose(fp_out);
  return 0;
}
