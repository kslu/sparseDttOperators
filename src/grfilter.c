#include "grfilter.h"
#include "matrices.h"

// #define SQRT2 1.4142135624
// cos(pi/4)
#define COSPI1_4 0.7071067812
// cos(k*pi/8) / sqrt(2), for k=1 and 3
#define COSPI1_8 0.6532814824
#define COSPI2_8 0.5000000000
#define COSPI3_8 0.2705980501
// cos(k*pi/16) / 2, for k=1, 3, 5, and 7
#define COSPI1_16 0.4903926402
#define COSPI3_16 0.4157348062
#define COSPI5_16 0.2777851165
#define COSPI7_16 0.0975451610

#if CONFIG_DEBUG
void show_coefficients(const double *input, const double *output, int n) {
  /*
  fprintf(stderr, "(");
  for (int i = 0; i < n - 1; i++)
    fprintf(stderr, "%.4lf,", input[i]);
  fprintf(stderr, "%.4lf)\n[", input[n - 1]);
  for (int i = 0; i < n - 1; i++)
    fprintf(stderr, "%.4lf,", output[i]);
  fprintf(stderr, "%.4lf]\n", output[n - 1]);
  */
}
#endif

void mat_times_vec(const double *input, double *output, const double *gftmtx,
                   int n) {
  for (int i = 0; i < n; i++) {
    output[i] = 0;
    for (int j = 0; j < n; j++)
      output[i] += input[j] * gftmtx[i * n + j];
  }
#if CONFIG_DEBUG
  show_coefficients(input, output, 10);
#endif
}

void dct4(const double *input, double *output) {
  double bf0[4] = {0};
  bf0[0] = input[0] + input[3];
  bf0[1] = input[1] + input[2];
  bf0[2] = input[1] - input[2];
  bf0[3] = input[0] - input[3];
  output[0] = COSPI2_8 * (bf0[0] + bf0[1]);
  output[1] = COSPI3_8 * bf0[2] + COSPI1_8 * bf0[3];
  output[2] = COSPI2_8 * (bf0[0] - bf0[1]);
  output[3] = -COSPI1_8 * bf0[2] + COSPI3_8 * bf0[3];
  return;
}

void idct4(const double *input, double *output) {
  double bf0[4] = {0};
  bf0[0] = COSPI2_8 * (input[0] + input[2]);
  bf0[1] = COSPI2_8 * (input[0] - input[2]);
  bf0[2] = COSPI3_8 * input[1] - COSPI1_8 * input[3];
  bf0[3] = COSPI1_8 * input[1] + COSPI3_8 * input[3];
  output[0] = bf0[0] + bf0[3];
  output[1] = bf0[1] + bf0[2];
  output[2] = bf0[1] - bf0[2];
  output[3] = bf0[0] - bf0[3];
  return;
}

void dct8(const double *input, double *output) {
  double bf0[8] = {0}, bf1[8] = {0};
  // stage 1
  bf0[0] = input[0] + input[7];
  bf0[1] = input[1] + input[6];
  bf0[2] = input[2] + input[5];
  bf0[3] = input[3] + input[4];
  bf0[4] = input[3] - input[4];
  bf0[5] = input[2] - input[5];
  bf0[6] = input[1] - input[6];
  bf0[7] = input[0] - input[7];
  // stage 2
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = bf0[1] - bf0[2];
  bf1[3] = bf0[0] - bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = COSPI1_4 * (-bf0[5] + bf0[6]);
  bf1[6] = COSPI1_4 * (bf0[5] + bf0[6]);
  bf1[7] = bf0[7];
  // stage 3
  bf0[0] = COSPI2_8 * (bf1[0] + bf1[1]);
  bf0[1] = COSPI2_8 * (bf1[0] - bf1[1]);
  bf0[2] = COSPI3_8 * bf1[2] + COSPI1_8 * bf1[3];
  bf0[3] = -COSPI1_8 * bf1[2] + COSPI3_8 * bf1[3];
  bf0[4] = bf1[4] + bf1[5];
  bf0[5] = bf1[4] - bf1[5];
  bf0[6] = -bf1[6] + bf1[7];
  bf0[7] = bf1[6] + bf1[7];
  // stage 4
  output[0] = bf0[0] * COSPI1_4;
  output[1] = COSPI7_16 * bf0[4] + COSPI1_16 * bf0[7];
  output[2] = bf0[2] * COSPI1_4;
  output[3] = -COSPI5_16 * bf0[5] + COSPI3_16 * bf0[6];
  output[4] = bf0[1] * COSPI1_4;
  output[5] = COSPI3_16 * bf0[5] + COSPI5_16 * bf0[6];
  output[6] = bf0[3] * COSPI1_4;
  output[7] = -COSPI1_16 * bf0[4] + COSPI7_16 * bf0[7];
  return;
}

void idct8(const double *input, double *output) {
  double bf0[8] = {0}, bf1[8] = {0};
  // stage 1
  bf0[0] = input[0] * COSPI1_4;
  bf0[1] = input[4] * COSPI1_4;
  bf0[2] = input[2] * COSPI1_4;
  bf0[3] = input[6] * COSPI1_4;
  bf0[4] = COSPI7_16 * input[1] - COSPI1_16 * input[7];
  bf0[5] = COSPI3_16 * input[5] - COSPI5_16 * input[3];
  bf0[6] = COSPI5_16 * input[5] + COSPI3_16 * input[3];
  bf0[7] = COSPI1_16 * input[1] + COSPI7_16 * input[7];
  // stage 2
  bf1[0] = COSPI2_8 * (bf0[0] + bf0[1]);
  bf1[1] = COSPI2_8 * (bf0[0] - bf0[1]);
  bf1[2] = COSPI3_8 * bf0[2] - COSPI1_8 * bf0[3];
  bf1[3] = COSPI1_8 * bf0[2] + COSPI3_8 * bf0[3];
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = bf0[4] - bf0[5];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[6] + bf0[7];
  // stage 3
  bf0[0] = bf1[0] + bf1[3];
  bf0[1] = bf1[1] + bf1[2];
  bf0[2] = bf1[1] - bf1[2];
  bf0[3] = bf1[0] - bf1[3];
  bf0[4] = bf1[4];
  bf0[5] = COSPI1_4 * (-bf1[5] + bf1[6]);
  bf0[6] = COSPI1_4 * (bf1[5] + bf1[6]);
  bf0[7] = bf1[7];
  // stage 4
  output[0] = bf0[0] + bf0[7];
  output[1] = bf0[1] + bf0[6];
  output[2] = bf0[2] + bf0[5];
  output[3] = bf0[3] + bf0[4];
  output[4] = bf0[3] - bf0[4];
  output[5] = bf0[2] - bf0[5];
  output[6] = bf0[1] - bf0[6];
  output[7] = bf0[0] - bf0[7];
  return;
}

void dct4x4(const double *input, double *output) {
  double bf0[4] = {0}, bf1[4] = {0};
  double temp[16] = {0};
  // apply dct4 row-wise
  for (int r = 0; r < 4; r++) {
    for (int i = 0; i < 4; i++)
      bf0[i] = input[4 * i + r];
    dct4(bf0, bf1);
    for (int i = 0; i < 4; i++)
      temp[4 * i + r] = bf1[i];
  }
  // apply dct4 column-wise
  for (int c = 0; c < 4; c++) {
    for (int i = 0; i < 4; i++)
      bf0[i] = temp[4 * c + i];
    dct4(bf0, bf1);
    for (int i = 0; i < 4; i++)
      output[4 * c + i] = bf1[i];
  }
}

void idct4x4(const double *input, double *output) {
  double bf0[4] = {0}, bf1[4] = {0};
  double temp[16] = {0};
  // apply idct4 row-wise
  for (int r = 0; r < 4; r++) {
    for (int i = 0; i < 4; i++)
      bf0[i] = input[4 * i + r];
    idct4(bf0, bf1);
    for (int i = 0; i < 4; i++)
      temp[4 * i + r] = bf1[i];
  }
  // apply idct4 column-wise
  for (int c = 0; c < 4; c++) {
    for (int i = 0; i < 4; i++)
      bf0[i] = temp[4 * c + i];
    idct4(bf0, bf1);
    for (int i = 0; i < 4; i++)
      output[4 * c + i] = bf1[i];
  }
}

void dct8x8(const double *input, double *output) {
  double bf0[8] = {0}, bf1[8] = {0};
  double temp[64] = {0};
  // apply dct8 row-wise
  for (int r = 0; r < 8; r++) {
    for (int i = 0; i < 8; i++)
      bf0[i] = input[8 * i + r];
    dct8(bf0, bf1);
    for (int i = 0; i < 8; i++)
      temp[8 * i + r] = bf1[i];
  }
  // apply dct8 column-wise
  for (int c = 0; c < 8; c++) {
    for (int i = 0; i < 8; i++)
      bf0[i] = temp[8 * c + i];
    dct8(bf0, bf1);
    for (int i = 0; i < 8; i++)
      output[8 * c + i] = bf1[i];
  }
}

void idct8x8(const double *input, double *output) {
  double bf0[8] = {0}, bf1[8] = {0};
  double temp[64] = {0};
  // apply idct8 row-wise
  for (int r = 0; r < 8; r++) {
    for (int i = 0; i < 8; i++)
      bf0[i] = input[8 * i + r];
    idct8(bf0, bf1);
    for (int i = 0; i < 8; i++)
      temp[8 * i + r] = bf1[i];
  }
  // apply idct8 column-wise
  for (int c = 0; c < 8; c++) {
    for (int i = 0; i < 8; i++)
      bf0[i] = temp[8 * c + i];
    idct8(bf0, bf1);
    for (int i = 0; i < 8; i++)
      output[8 * c + i] = bf1[i];
  }
}

void apply_sparse_laplacian(const double *input, double *output, int n,
                            const int nedges, const double mev,
                            const int *adjlist, const double *wlist) {
  int s = 0, t = 0;
  for (int i = 0; i < n; i++)
    output[i] = -mev * input[i];
  for (int i = 0; i < nedges; i++) {
    s = adjlist[i * 2];
    t = adjlist[i * 2 + 1];
    if (s == t) {
      output[s] += wlist[i] * input[s];
    } else {
      output[s] += wlist[i] * (input[s] - input[t]);
      output[t] += wlist[i] * (input[t] - input[s]);
    }
  }
  return;
}

void apply_sparse_operator(const double *input, double *output, int n,
                           const int nedges, const double mev,
                           const int *adjlist, const double *wlist) {
  int s = 0, t = 0;
  for (int i = 0; i < n; i++)
    output[i] = -mev * input[i];
  for (int i = 0; i < nedges; i++) {
    s = adjlist[i * 2];
    t = adjlist[i * 2 + 1];
    if (s == t) {
      output[s] += wlist[i] * input[t];
    } else {
      output[s] += wlist[i] * input[t];
      output[t] += wlist[i] * input[s];
    }
  }
  return;
}

void fir_graph_filter(const double *input, double *output, int n, int order,
                      const double *coeffs, const int nedges, const double mev,
                      const int *adjlist, const double *wlist) {
  double temp[MAX_GRAPH_SIZE];
  for (int i = 0; i < n; i++)
    output[i] = coeffs[order] * input[i];
  for (int i = order - 1; i >= 0; i--) {
    apply_sparse_laplacian(output, temp, n, nedges, mev, adjlist, wlist);
    for (int j = 0; j < n; j++)
      output[j] = temp[j] + coeffs[i] * input[j];
  }
  return;
}

void get_multishift_terms(const int *powers, int ord, int m, int nops,
                          int *idx_list, int *pow_list) {
  for (int k = 0; k < m; k++) {
    int l = 0;
    for (int opid = 0; opid < nops; opid++) {
      int op_pow = powers[k * nops + opid];
      if (op_pow == 0)
        continue;
      idx_list[k * ord + l] = opid;
      pow_list[k * ord + l] = op_pow;
      l++;
    }
  }
  return;
}

void multishift_graph_filter(const double *input, double *output, int n,
                             int ord, int m, const double *coeffs,
                             const int *idx_list, const int *pow_list,
                             const int *nedges, const int *alists[],
                             const double *wlists[]) {
  double temp1[MAX_GRAPH_SIZE], temp2[MAX_GRAPH_SIZE];
  for (int i = 0; i < n; i++)
    output[i] = 0;
  for (int k = 0; k < m; k++) {
    for (int j = 0; j < n; j++)
      temp2[j] = input[j];
    for (int op = 0; op < ord; op++) {
      int pow = pow_list[k * ord + op], idx = idx_list[k * ord + op];
      if (pow == 0)
        break;
      if (idx == 0)
        continue;
      for (int j = 0; j < pow; j++) {
        apply_sparse_operator(temp2, temp1, n, nedges[idx], 0, alists[idx],
                              wlists[idx]);
        for (int l = 0; l < n; l++)
          temp2[l] = temp1[l];
      }
    }
    for (int j = 0; j < n; j++)
      output[j] += coeffs[k] * temp2[j];
  }
  return;
}
