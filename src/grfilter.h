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

// cospi[k] = cos(k*pi/128)
static const double cospi[64] = {
    1.0000000000, 0.9996988187, 0.9987954562, 0.9972904567, 0.9951847267,
    0.9924795346, 0.9891765100, 0.9852776424, 0.9807852804, 0.9757021300,
    0.9700312532, 0.9637760658, 0.9569403357, 0.9495281806, 0.9415440652,
    0.9329927988, 0.9238795325, 0.9142097557, 0.9039892931, 0.8932243012,
    0.8819212643, 0.8700869911, 0.8577286100, 0.8448535652, 0.8314696123,
    0.8175848132, 0.8032075315, 0.7883464276, 0.7730104534, 0.7572088465,
    0.7409511254, 0.7242470830, 0.7071067812, 0.6895405447, 0.6715589548,
    0.6531728430, 0.6343932842, 0.6152315906, 0.5956993045, 0.5758081914,
    0.5555702330, 0.5349976199, 0.5141027442, 0.4928981922, 0.4713967368,
    0.4496113297, 0.4275550934, 0.4052413140, 0.3826834324, 0.3598950365,
    0.3368898534, 0.3136817404, 0.2902846773, 0.2667127575, 0.2429801799,
    0.2191012402, 0.1950903220, 0.1709618888, 0.1467304745, 0.1224106752,
    0.0980171403, 0.0735645636, 0.0490676743, 0.0245412285,
};

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

static const double h32_lp[32] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
};

static const double h4x4_tik[16] = {
    1.0000000000, 0.9646814192, 0.8888888889, 0.8241384565,
    0.9646814192, 0.9317725357, 0.8608729070, 0.8000000000,
    0.8888888889, 0.8608729070, 0.8000000000, 0.7471672940,
    0.8241384565, 0.8000000000, 0.7471672940, 0.7008805255,
};

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

static const double h32_tik[32] = {
    1.0000000000, 0.9993984529, 0.9976039151, 0.9946463577, 0.9905746246,
    0.9854548423, 0.9793683344, 0.9724091607, 0.9646814192, 0.9562964495,
    0.9473700677, 0.9380199492, 0.9283632483, 0.9185145213, 0.9085839938,
    0.8986761855, 0.8888888889, 0.8793124784, 0.8700295179, 0.8611146244,
    0.8526345430, 0.8446483895, 0.8372080164, 0.8303584631, 0.8241384565,
    0.8185809315, 0.8137135459, 0.8095591724, 0.8061363476, 0.8034596704,
    0.8015401369, 0.8003854075,
};

void exact_filter_4x4(const double *input, double *output, const double *h);
void exact_filter_8x8(const double *input, double *output, const double *h);
void exact_filter_32(const double *input, double *output, const double *h);
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
