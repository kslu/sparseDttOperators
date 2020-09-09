#include "grfilter.h"
#include "dct.h"
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

void exact_filter_4x4(const double *input, double *output, const double *h) {
  double temp[16] = {0};
  dct4x4(input, temp);
  for (int i = 0; i < 16; i++)
    temp[i] *= h[i];
  idct4x4(temp, output);
  return;
}

void exact_filter_8x8(const double *input, double *output, const double *h) {
  double temp[64] = {0};
  dct8x8(input, temp);
  for (int i = 0; i < 64; i++)
    temp[i] *= h[i];
  idct8x8(temp, output);
  return;
}

void exact_filter_16x16(const double *input, double *output, const double *h) {
  double temp[256] = {0};
  dct16x16(input, temp);
  for (int i = 0; i < 256; i++)
    temp[i] *= h[i];
  idct16x16(temp, output);
  return;
}

void exact_filter_32(const double *input, double *output, const double *h) {
  double temp[32] = {0};
  dct32(input, temp);
  for (int i = 0; i < 32; i++)
    temp[i] *= h[i];
  idct32(temp, output);
  return;
}

void exact_filter_64(const double *input, double *output, const double *h) {
  double temp[64] = {0};
  dct64(input, temp);
  for (int i = 0; i < 64; i++)
    temp[i] *= h[i];
  idct64(temp, output);
  return;
}

// y = (L - mev * I) * x
void apply_sparse_laplacian(const double *input, double *output, int n,
                            const int nedges, const double mev,
                            const int *adjlist, const double *wlist) {
  int s = 0, t = 0;
  if (mev == 0) {
    memset(output, 0, n * sizeof(double));
  } else {
    for (int i = 0; i < n; i++)
      output[i] = -mev * input[i];
  }
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
  if (mev == 0) {
    memset(output, 0, n * sizeof(double));
  } else {
    for (int i = 0; i < n; i++)
      output[i] = -mev * input[i];
  }
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

void buffer_add(const double *input, double *output, int dim, int len) {
  int ind=0;
  for (int j=0; j<dim; j++)
    for (int i=0; i<len; i++)
      output[i] += input[ind++];
  return;
}

// Polynomial graph filter (polynomial of L)
// y = (c[0] + c[1] * L + c[2] * L^2 + ... + c[K-1] * L^{K-1}) * x
void pgf(const double *input, double *output, int n, int order,
         const double *coeffs, const int nedges, const double mev,
         const int *adjlist, const double *wlist) {
  double temp[MAX_GRAPH_SIZE];
  for (int i = 0; i < n; i++)
    output[i] = coeffs[order] * input[i];
  for (int i = order - 1; i >= 0; i--) {
    apply_sparse_laplacian(output, temp, n, nedges, mev, adjlist, wlist);
    memcpy(output, temp, n * sizeof(double));
    if (coeffs[i] != 0) {
      for (int j = 0; j < n; j++)
        output[j] += coeffs[i] * input[j];
    }
  }
  return;
}

// Polynomial graph filter (polynomial of S)
// y = (c[0] + c[1] * S + c[2] * S^2 + ... + c[K-1] * S^{K-1}) * x
void pgf_s(const double *input, double *output, int n, int order,
           const double *coeffs, const int nedges, const double mev,
           const int *adjlist, const double *wlist) {
  double temp[MAX_GRAPH_SIZE];
  for (int i = 0; i < n; i++)
    output[i] = coeffs[order] * input[i];
  for (int i = order - 1; i >= 0; i--) {
    apply_sparse_operator(output, temp, n, nedges, mev, adjlist, wlist);
    memcpy(output, temp, n * sizeof(double));
    if (coeffs[i] != 0) {
      for (int j = 0; j < n; j++)
        output[j] += coeffs[i] * input[j];
    }
  }
  return;
}

// Chebyshev graph filter
void chebyshev_gf(const double *input, double *output, int n, int order,
                  const double *coeffs, const int nedges, const double mev,
                  const int *adjlist, const double *wlist) {
  // twf_new: T_k(L) f
  // twf_cur: T_{k-1}(L) f
  // twf_old: T_{k-2}(L) f
  double twf_old[MAX_GRAPH_SIZE], twf_cur[MAX_GRAPH_SIZE],
      twf_new[MAX_GRAPH_SIZE];
  double inv_mev = 1.0 / mev;
  memcpy(twf_old, input, n * sizeof(double));
  apply_sparse_laplacian(input, twf_cur, n, nedges, mev, adjlist, wlist);

  for (int i = 0; i < n; i++) {
    twf_cur[i] = inv_mev * twf_cur[i];
    output[i] = 0.5 * coeffs[0] * twf_old[i] + coeffs[1] * twf_cur[i];
  }

  // TODO: check "k<=order-1" or "k<=order"
  for (int k = 2; k <= order - 1; k++) {
    // twf_new is used as a buffer here
    apply_sparse_laplacian(twf_cur, twf_new, n, nedges, mev, adjlist, wlist);
    for (int i = 0; i < n; i++) {
      twf_new[i] = 2 * inv_mev * twf_new[i] - twf_old[i];
      output[i] += coeffs[k] * twf_new[i];
    }
    if (k < order - 1) {
      memcpy(twf_old, twf_cur, n * sizeof(double));
      memcpy(twf_cur, twf_new, n * sizeof(double));
    }
  }

  return;
}

// ARMA graph filter, conjugate gradient (CG) implementation
// y = Hx, where the frequency response of H has the form
//         b[0] + b[1]*z + ... + b[kb-1]*z^{kb-1}
// h(z) = ----------------------------------------
//         a[0] + a[1]*z + ... + a[ka-1]*z^{ka-1}
// with coefficients a[i], b[i] chosen to approximate the desired frequency
// response
void armagf_cg(const double *input, double *output, int n, int tmax, int kb,
               const double *b, int ka, const double *a, const int nedges,
               const int *adjlist, const double *wlist) {
  double lp[MAX_GRAPH_SIZE], ay[MAX_GRAPH_SIZE], r[MAX_GRAPH_SIZE],
      p[MAX_GRAPH_SIZE];
  double rsold = 0, dot_temp = 0, alpha = 0;
  pgf(input, output, n, kb, b, nedges, 0, adjlist, wlist); // y = b(L) * x
  pgf(output, ay, n, ka, a, nedges, 0, adjlist, wlist);    // ay = a(L) * y

  for (int i = 0; i < n; i++) {
    r[i] = output[i] - ay[i]; // r = b(L) * x - a(L) * y
    rsold += r[i] * r[i];
  }
  memcpy(p, r, n * sizeof(double));

  for (int k = 0; k < tmax; k++) {
    pgf(p, lp, n, ka, a, nedges, 0, adjlist, wlist); // lp = a(L) * p

    dot_temp = 0;
    for (int i = 0; i < n; i++)
      dot_temp += p[i] * lp[i];
    alpha = rsold / dot_temp;

    dot_temp = 0;
    for (int i = 0; i < n; i++) {
      output[i] += alpha * p[i];
      r[i] -= alpha * lp[i];
      dot_temp += r[i] * r[i];
    }

    if (dot_temp <= 1e-12)
      break;

    for (int i = 0; i < n; i++)
      p[i] = r[i] + (dot_temp / rsold) * p[i];
    rsold = dot_temp;
  }

  return;
}

void get_mpgf_terms(const int *powers, int ord, int m, int nops, int *idx_list,
                    int *pow_list) {
  int l = 0, op_pow = 0;
  for (int k = 0; k < m; k++) {
    l = 0;
    for (int opid = 0; opid < nops; opid++) {
      op_pow = powers[k * nops + opid];
      if (op_pow == 0)
        continue;
      idx_list[k * ord + l] = opid;
      pow_list[k * ord + l] = op_pow;
      l++;
    }
  }
  return;
}

// Multivariate polynomial graph filter
void mpgf(const double *input, double *output, int n, int ord, int m,
          const double *coeffs, const int *idx_list, const int *pow_list,
          const int *nedges, const int *alists[], const double *wlists[]) {
  double temp1[MAX_GRAPH_SIZE] = {0}, temp2[MAX_GRAPH_SIZE] = {0};
  int pow = 0, idx = 0, idxlist = 0;
  for (int k = 0; k < m; k++) {
    for (int j = 0; j < n; j++)
      temp2[j] = input[j];
    for (int op = 0; op < ord; op++) {
      pow = pow_list[idxlist + op];
      idx = idx_list[idxlist + op];
      if (pow == 0)
        break;
      if (idx == 0)
        continue;
      for (int j = 0; j < pow; j++) {
        apply_sparse_operator(temp2, temp1, n, nedges[idx], 0, alists[idx],
                              wlists[idx]);
        memcpy(temp2, temp1, n * sizeof(double));
      }
    }
    idxlist += ord;
    if (coeffs[k] != 0) {
      for (int j = 0; j < n; j++)
        output[j] += coeffs[k] * temp2[j];
    }
  }
  return;
}
