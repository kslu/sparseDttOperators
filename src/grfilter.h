#ifndef GRFILTER_H_
#define GRFILTER_H_

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define CONFIG_DEBUG 0
#define MAX_NUM_DATA 2000
#define MAX_GRAPH_SIZE 128
#define MAX_EDGE_NUMBER MAX_GRAPH_SIZE *MAX_GRAPH_SIZE
#define BATCH_SIZE 500

#define SQRT2 1.414213564
#define INVSQRT2 0.7071067812

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
void multishift_graph_filter(const double *input, double *output, int n, int m,
                             const double *coeffs, int nops, const int *nedges,
                             const int *powers, const int *alists[],
                             const double *wlists[]);
#endif
