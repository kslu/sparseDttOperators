#include "grfilter.h"
#include "matrices.h"

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
