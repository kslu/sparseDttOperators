#ifndef GRFILTER_H_
#define GRFILTER_H_

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define CONFIG_DEBUG 0
#define MAX_GRAPH_SIZE 128
#define MAX_EDGE_NUMBER MAX_GRAPH_SIZE *MAX_GRAPH_SIZE
#define BATCH_SIZE 100

#if CONFIG_DEBUG
void show_coefficients(const double *input, const double *output, int n);
#endif

void mat_times_vec(const double *input, double *output, const double *gftmtx,
                   int n);
void apply_sparse_operator(const double *input, double *output, int n,
                           const int nedges, const double mev,
                           const int *adjlist, const double *wlist);
void fir_graph_filter(const double *input, double *output, int n, int order,
                      const double *coeffs, const int nedges, const double mev,
                      const int *adjlist, const double *wlist);
void get_multishift_terms(const int *powers, int ord, int m, int nops,
                          int *idx_list, int *pow_list);
void multishift_graph_filter(const double *input, double *output, int n,
                             int ord, int m, const double *coeffs,
                             const int *idx_list, const int *pow_list,
                             const int *nedges, const int *alists[],
                             const double *wlists[]);
#endif
