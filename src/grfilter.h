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

void dct4(const double *input, double *output);
void idct4(const double *input, double *output);
void dct8(const double *input, double *output);
void idct8(const double *input, double *output);
void dct4x4(const double *input, double *output);
void idct4x4(const double *input, double *output);
void dct8x8(const double *input, double *output);
void idct8x8(const double *input, double *output);

static const double h4x4_lp[16] = {
    1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};

static const double h8x8_lp[64] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};

static const double h4x4_tik[16] = {
    1.0000000000, 0.9646814192, 0.8888888889, 0.8241384565,
    0.9646814192, 0.9317725357, 0.8608729070, 0.8000000000,
    0.8888888889, 0.8608729070, 0.8000000000, 0.7471672940,
    0.8241384565, 0.8000000000, 0.7471672940, 0.7008805255};

static const double h8x8_tik[64] = {
    1.0000000000, 0.9905746246, 0.9646814192, 0.9283632483, 0.8888888889,
    0.8526345430, 0.8241384565, 0.8061363476, 0.9905746246, 0.9813252655,
    0.9559071475, 0.9202344214, 0.8814338713, 0.8457728995, 0.8177260731,
    0.8000000000, 0.9646814192, 0.9559071475, 0.9317725357, 0.8978464247,
    0.8608729070, 0.8268241045, 0.8000000000, 0.7830261272, 0.9283632483,
    0.9202344214, 0.8978464247, 0.8663040408, 0.8318328656, 0.8000000000,
    0.7748616780, 0.7589271582, 0.8888888889, 0.8814338713, 0.8608729070,
    0.8318328656, 0.8000000000, 0.7705137166, 0.7471672940, 0.7323405550,
    0.8526345430, 0.8457728995, 0.8268241045, 0.8000000000, 0.7705137166,
    0.7431237691, 0.7213842504, 0.7075536593, 0.8241384565, 0.8177260731,
    0.8000000000, 0.7748616780, 0.7471672940, 0.7213842504, 0.7008805255,
    0.6878178500, 0.8061363476, 0.8000000000, 0.7830261272, 0.7589271582,
    0.7323405550, 0.7075536593, 0.6878178500, 0.6752331775,
};

void exact_filter_4x4(const double *input, double *output, double *h);
void exact_filter_8x8(const double *input, double *output, double *h);
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
