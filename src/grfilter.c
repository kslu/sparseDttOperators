#include "grfilter.h"
#include "matrices.h"

#define SQRT2 1.4142135624
#define INVSQRT2 0.7071067812

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
  output[0] = cospi[32] * (bf0[0] + bf0[1]) * INVSQRT2;
  output[1] = (cospi[48] * bf0[2] + cospi[16] * bf0[3]) * INVSQRT2;
  output[2] = cospi[32] * (bf0[0] - bf0[1]) * INVSQRT2;
  output[3] = (-cospi[16] * bf0[2] + cospi[48] * bf0[3]) * INVSQRT2;
  return;
}

void idct4(const double *input, double *output) {
  double bf0[4] = {0};
  bf0[0] = cospi[32] * (input[0] + input[2]) * INVSQRT2;
  bf0[1] = cospi[32] * (input[0] - input[2]) * INVSQRT2;
  bf0[2] = (cospi[48] * input[1] - cospi[16] * input[3]) * INVSQRT2;
  bf0[3] = (cospi[16] * input[1] + cospi[48] * input[3]) * INVSQRT2;
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
  bf1[5] = cospi[32] * (-bf0[5] + bf0[6]);
  bf1[6] = cospi[32] * (bf0[5] + bf0[6]);
  bf1[7] = bf0[7];
  // stage 3
  bf0[0] = cospi[32] * (bf1[0] + bf1[1]);
  bf0[1] = cospi[32] * (bf1[0] - bf1[1]);
  bf0[2] = cospi[48] * bf1[2] + cospi[16] * bf1[3];
  bf0[3] = -cospi[16] * bf1[2] + cospi[48] * bf1[3];
  bf0[4] = bf1[4] + bf1[5];
  bf0[5] = bf1[4] - bf1[5];
  bf0[6] = -bf1[6] + bf1[7];
  bf0[7] = bf1[6] + bf1[7];
  // stage 4
  output[0] = bf0[0] * cospi[32] * INVSQRT2;
  output[1] = (cospi[56] * bf0[4] + cospi[8] * bf0[7]) / 2;
  output[2] = bf0[2] * cospi[32] * INVSQRT2;
  output[3] = (-cospi[40] * bf0[5] + cospi[24] * bf0[6]) / 2;
  output[4] = bf0[1] * cospi[32] * INVSQRT2;
  output[5] = (cospi[24] * bf0[5] + cospi[40] * bf0[6]) / 2;
  output[6] = bf0[3] * cospi[32] * INVSQRT2;
  output[7] = (-cospi[8] * bf0[4] + cospi[56] * bf0[7]) / 2;
  return;
}

void idct8(const double *input, double *output) {
  double bf0[8] = {0}, bf1[8] = {0};
  // stage 1
  bf0[0] = input[0] * cospi[32] * INVSQRT2;
  bf0[1] = input[4] * cospi[32] * INVSQRT2;
  bf0[2] = input[2] * cospi[32] * INVSQRT2;
  bf0[3] = input[6] * cospi[32] * INVSQRT2;
  bf0[4] = (cospi[56] * input[1] - cospi[8] * input[7]) / 2;
  bf0[5] = (cospi[24] * input[5] - cospi[40] * input[3]) / 2;
  bf0[6] = (cospi[40] * input[5] + cospi[24] * input[3]) / 2;
  bf0[7] = (cospi[8] * input[1] + cospi[56] * input[7]) / 2;
  // stage 2
  bf1[0] = cospi[32] * (bf0[0] + bf0[1]);
  bf1[1] = cospi[32] * (bf0[0] - bf0[1]);
  bf1[2] = cospi[48] * bf0[2] - cospi[16] * bf0[3];
  bf1[3] = cospi[16] * bf0[2] + cospi[48] * bf0[3];
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
  bf0[5] = cospi[32] * (-bf1[5] + bf1[6]);
  bf0[6] = cospi[32] * (bf1[5] + bf1[6]);
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

void dct32(const double *input, double *output) {
  double bf0[32] = {0}, bf1[32] = {0};
  // stage 1;
  bf0[0] = input[0] + input[31];
  bf0[1] = input[1] + input[30];
  bf0[2] = input[2] + input[29];
  bf0[3] = input[3] + input[28];
  bf0[4] = input[4] + input[27];
  bf0[5] = input[5] + input[26];
  bf0[6] = input[6] + input[25];
  bf0[7] = input[7] + input[24];
  bf0[8] = input[8] + input[23];
  bf0[9] = input[9] + input[22];
  bf0[10] = input[10] + input[21];
  bf0[11] = input[11] + input[20];
  bf0[12] = input[12] + input[19];
  bf0[13] = input[13] + input[18];
  bf0[14] = input[14] + input[17];
  bf0[15] = input[15] + input[16];
  bf0[16] = -input[16] + input[15];
  bf0[17] = -input[17] + input[14];
  bf0[18] = -input[18] + input[13];
  bf0[19] = -input[19] + input[12];
  bf0[20] = -input[20] + input[11];
  bf0[21] = -input[21] + input[10];
  bf0[22] = -input[22] + input[9];
  bf0[23] = -input[23] + input[8];
  bf0[24] = -input[24] + input[7];
  bf0[25] = -input[25] + input[6];
  bf0[26] = -input[26] + input[5];
  bf0[27] = -input[27] + input[4];
  bf0[28] = -input[28] + input[3];
  bf0[29] = -input[29] + input[2];
  bf0[30] = -input[30] + input[1];
  bf0[31] = -input[31] + input[0];

  // stage 2
  bf1[0] = bf0[0] + bf0[15];
  bf1[1] = bf0[1] + bf0[14];
  bf1[2] = bf0[2] + bf0[13];
  bf1[3] = bf0[3] + bf0[12];
  bf1[4] = bf0[4] + bf0[11];
  bf1[5] = bf0[5] + bf0[10];
  bf1[6] = bf0[6] + bf0[9];
  bf1[7] = bf0[7] + bf0[8];
  bf1[8] = -bf0[8] + bf0[7];
  bf1[9] = -bf0[9] + bf0[6];
  bf1[10] = -bf0[10] + bf0[5];
  bf1[11] = -bf0[11] + bf0[4];
  bf1[12] = -bf0[12] + bf0[3];
  bf1[13] = -bf0[13] + bf0[2];
  bf1[14] = -bf0[14] + bf0[1];
  bf1[15] = -bf0[15] + bf0[0];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = cospi[32] * (-bf0[20] + bf0[27]);
  bf1[21] = cospi[32] * (-bf0[21] + bf0[26]);
  bf1[22] = cospi[32] * (-bf0[22] + bf0[25]);
  bf1[23] = cospi[32] * (-bf0[23] + bf0[24]);
  bf1[24] = cospi[32] * (bf0[24] + bf0[23]);
  bf1[25] = cospi[32] * (bf0[25] + bf0[22]);
  bf1[26] = cospi[32] * (bf0[26] + bf0[21]);
  bf1[27] = cospi[32] * (bf0[27] + bf0[20]);
  bf1[28] = bf0[28];
  bf1[29] = bf0[29];
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 3
  bf0[0] = bf1[0] + bf1[7];
  bf0[1] = bf1[1] + bf1[6];
  bf0[2] = bf1[2] + bf1[5];
  bf0[3] = bf1[3] + bf1[4];
  bf0[4] = -bf1[4] + bf1[3];
  bf0[5] = -bf1[5] + bf1[2];
  bf0[6] = -bf1[6] + bf1[1];
  bf0[7] = -bf1[7] + bf1[0];
  bf0[8] = bf1[8];
  bf0[9] = bf1[9];
  bf0[10] = cospi[32] * (-bf1[10] + bf1[13]);
  bf0[11] = cospi[32] * (-bf1[11] + bf1[12]);
  bf0[12] = cospi[32] * (bf1[12] + bf1[11]);
  bf0[13] = cospi[32] * (bf1[13] + bf1[10]);
  bf0[14] = bf1[14];
  bf0[15] = bf1[15];
  bf0[16] = bf1[16] + bf1[23];
  bf0[17] = bf1[17] + bf1[22];
  bf0[18] = bf1[18] + bf1[21];
  bf0[19] = bf1[19] + bf1[20];
  bf0[20] = -bf1[20] + bf1[19];
  bf0[21] = -bf1[21] + bf1[18];
  bf0[22] = -bf1[22] + bf1[17];
  bf0[23] = -bf1[23] + bf1[16];
  bf0[24] = -bf1[24] + bf1[31];
  bf0[25] = -bf1[25] + bf1[30];
  bf0[26] = -bf1[26] + bf1[29];
  bf0[27] = -bf1[27] + bf1[28];
  bf0[28] = bf1[28] + bf1[27];
  bf0[29] = bf1[29] + bf1[26];
  bf0[30] = bf1[30] + bf1[25];
  bf0[31] = bf1[31] + bf1[24];

  // stage 4
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = -bf0[2] + bf0[1];
  bf1[3] = -bf0[3] + bf0[0];
  bf1[4] = bf0[4];
  bf1[5] = cospi[32] * (-bf0[5] + bf0[6]);
  bf1[6] = cospi[32] * (bf0[6] + bf0[5]);
  bf1[7] = bf0[7];
  bf1[8] = bf0[8] + bf0[11];
  bf1[9] = bf0[9] + bf0[10];
  bf1[10] = -bf0[10] + bf0[9];
  bf1[11] = -bf0[11] + bf0[8];
  bf1[12] = -bf0[12] + bf0[15];
  bf1[13] = -bf0[13] + bf0[14];
  bf1[14] = bf0[14] + bf0[13];
  bf1[15] = bf0[15] + bf0[12];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = -cospi[16] * bf0[18] + cospi[48] * bf0[29];
  bf1[19] = -cospi[16] * bf0[19] + cospi[48] * bf0[28];
  bf1[20] = -cospi[48] * bf0[20] - cospi[16] * bf0[27];
  bf1[21] = -cospi[48] * bf0[21] - cospi[16] * bf0[26];
  bf1[22] = bf0[22];
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = cospi[48] * bf0[26] - cospi[16] * bf0[21];
  bf1[27] = cospi[48] * bf0[27] - cospi[16] * bf0[20];
  bf1[28] = cospi[16] * bf0[28] + cospi[48] * bf0[19];
  bf1[29] = cospi[16] * bf0[29] + cospi[48] * bf0[18];
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 5
  bf0[0] = cospi[32] * (bf1[0] + bf1[1]);
  bf0[1] = cospi[32] * (-bf1[1] + bf1[0]);
  bf0[2] = cospi[48] * bf1[2] + cospi[16] * bf1[3];
  bf0[3] = cospi[48] * bf1[3] - cospi[16] * bf1[2];
  bf0[4] = bf1[4] + bf1[5];
  bf0[5] = -bf1[5] + bf1[4];
  bf0[6] = -bf1[6] + bf1[7];
  bf0[7] = bf1[7] + bf1[6];
  bf0[8] = bf1[8];
  bf0[9] = -cospi[16] * bf1[9] + cospi[48] * bf1[14];
  bf0[10] = -cospi[48] * bf1[10] - cospi[16] * bf1[13];
  bf0[11] = bf1[11];
  bf0[12] = bf1[12];
  bf0[13] = cospi[48] * bf1[13] - cospi[16] * bf1[10];
  bf0[14] = cospi[16] * bf1[14] + cospi[48] * bf1[9];
  bf0[15] = bf1[15];
  bf0[16] = bf1[16] + bf1[19];
  bf0[17] = bf1[17] + bf1[18];
  bf0[18] = -bf1[18] + bf1[17];
  bf0[19] = -bf1[19] + bf1[16];
  bf0[20] = -bf1[20] + bf1[23];
  bf0[21] = -bf1[21] + bf1[22];
  bf0[22] = bf1[22] + bf1[21];
  bf0[23] = bf1[23] + bf1[20];
  bf0[24] = bf1[24] + bf1[27];
  bf0[25] = bf1[25] + bf1[26];
  bf0[26] = -bf1[26] + bf1[25];
  bf0[27] = -bf1[27] + bf1[24];
  bf0[28] = -bf1[28] + bf1[31];
  bf0[29] = -bf1[29] + bf1[30];
  bf0[30] = bf1[30] + bf1[29];
  bf0[31] = bf1[31] + bf1[28];

  // stage 6
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = cospi[56] * bf0[4] + cospi[8] * bf0[7];
  bf1[5] = cospi[24] * bf0[5] + cospi[40] * bf0[6];
  bf1[6] = cospi[24] * bf0[6] - cospi[40] * bf0[5];
  bf1[7] = cospi[56] * bf0[7] - cospi[8] * bf0[4];
  bf1[8] = bf0[8] + bf0[9];
  bf1[9] = -bf0[9] + bf0[8];
  bf1[10] = -bf0[10] + bf0[11];
  bf1[11] = bf0[11] + bf0[10];
  bf1[12] = bf0[12] + bf0[13];
  bf1[13] = -bf0[13] + bf0[12];
  bf1[14] = -bf0[14] + bf0[15];
  bf1[15] = bf0[15] + bf0[14];
  bf1[16] = bf0[16];
  bf1[17] = -cospi[8] * bf0[17] + cospi[56] * bf0[30];
  bf1[18] = -cospi[56] * bf0[18] - cospi[8] * bf0[29];
  bf1[19] = bf0[19];
  bf1[20] = bf0[20];
  bf1[21] = -cospi[40] * bf0[21] + cospi[24] * bf0[26];
  bf1[22] = -cospi[24] * bf0[22] - cospi[40] * bf0[25];
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = cospi[24] * bf0[25] - cospi[40] * bf0[22];
  bf1[26] = cospi[40] * bf0[26] + cospi[24] * bf0[21];
  bf1[27] = bf0[27];
  bf1[28] = bf0[28];
  bf1[29] = cospi[56] * bf0[29] - cospi[8] * bf0[18];
  bf1[30] = cospi[8] * bf0[30] + cospi[56] * bf0[17];
  bf1[31] = bf0[31];

  // stage 7
  bf0[0] = bf1[0];
  bf0[1] = bf1[1];
  bf0[2] = bf1[2];
  bf0[3] = bf1[3];
  bf0[4] = bf1[4];
  bf0[5] = bf1[5];
  bf0[6] = bf1[6];
  bf0[7] = bf1[7];
  bf0[8] = cospi[60] * bf1[8] + cospi[4] * bf1[15];
  bf0[9] = cospi[28] * bf1[9] + cospi[36] * bf1[14];
  bf0[10] = cospi[44] * bf1[10] + cospi[20] * bf1[13];
  bf0[11] = cospi[12] * bf1[11] + cospi[52] * bf1[12];
  bf0[12] = cospi[12] * bf1[12] - cospi[52] * bf1[11];
  bf0[13] = cospi[44] * bf1[13] - cospi[20] * bf1[10];
  bf0[14] = cospi[28] * bf1[14] - cospi[36] * bf1[9];
  bf0[15] = cospi[60] * bf1[15] - cospi[4] * bf1[8];
  bf0[16] = bf1[16] + bf1[17];
  bf0[17] = -bf1[17] + bf1[16];
  bf0[18] = -bf1[18] + bf1[19];
  bf0[19] = bf1[19] + bf1[18];
  bf0[20] = bf1[20] + bf1[21];
  bf0[21] = -bf1[21] + bf1[20];
  bf0[22] = -bf1[22] + bf1[23];
  bf0[23] = bf1[23] + bf1[22];
  bf0[24] = bf1[24] + bf1[25];
  bf0[25] = -bf1[25] + bf1[24];
  bf0[26] = -bf1[26] + bf1[27];
  bf0[27] = bf1[27] + bf1[26];
  bf0[28] = bf1[28] + bf1[29];
  bf0[29] = -bf1[29] + bf1[28];
  bf0[30] = -bf1[30] + bf1[31];
  bf0[31] = bf1[31] + bf1[30];

  // stage 8
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = cospi[62] * bf0[16] + cospi[2] * bf0[31];
  bf1[17] = cospi[30] * bf0[17] + cospi[34] * bf0[30];
  bf1[18] = cospi[46] * bf0[18] + cospi[18] * bf0[29];
  bf1[19] = cospi[14] * bf0[19] + cospi[50] * bf0[28];
  bf1[20] = cospi[54] * bf0[20] + cospi[10] * bf0[27];
  bf1[21] = cospi[22] * bf0[21] + cospi[42] * bf0[26];
  bf1[22] = cospi[38] * bf0[22] + cospi[26] * bf0[25];
  bf1[23] = cospi[6] * bf0[23] + cospi[58] * bf0[24];
  bf1[24] = cospi[6] * bf0[24] - cospi[58] * bf0[23];
  bf1[25] = cospi[38] * bf0[25] - cospi[26] * bf0[22];
  bf1[26] = cospi[22] * bf0[26] - cospi[42] * bf0[21];
  bf1[27] = cospi[54] * bf0[27] - cospi[10] * bf0[20];
  bf1[28] = cospi[14] * bf0[28] - cospi[50] * bf0[19];
  bf1[29] = cospi[46] * bf0[29] - cospi[18] * bf0[18];
  bf1[30] = cospi[30] * bf0[30] - cospi[34] * bf0[17];
  bf1[31] = cospi[62] * bf0[31] - cospi[2] * bf0[16];

  // stage 9
  output[0] = bf1[0] / 4;
  output[1] = bf1[16] / 4;
  output[2] = bf1[8] / 4;
  output[3] = bf1[24] / 4;
  output[4] = bf1[4] / 4;
  output[5] = bf1[20] / 4;
  output[6] = bf1[12] / 4;
  output[7] = bf1[28] / 4;
  output[8] = bf1[2] / 4;
  output[9] = bf1[18] / 4;
  output[10] = bf1[10] / 4;
  output[11] = bf1[26] / 4;
  output[12] = bf1[6] / 4;
  output[13] = bf1[22] / 4;
  output[14] = bf1[14] / 4;
  output[15] = bf1[30] / 4;
  output[16] = bf1[1] / 4;
  output[17] = bf1[17] / 4;
  output[18] = bf1[9] / 4;
  output[19] = bf1[25] / 4;
  output[20] = bf1[5] / 4;
  output[21] = bf1[21] / 4;
  output[22] = bf1[13] / 4;
  output[23] = bf1[29] / 4;
  output[24] = bf1[3] / 4;
  output[25] = bf1[19] / 4;
  output[26] = bf1[11] / 4;
  output[27] = bf1[27] / 4;
  output[28] = bf1[7] / 4;
  output[29] = bf1[23] / 4;
  output[30] = bf1[15] / 4;
  output[31] = bf1[31] / 4;
  return;
}

void idct32(const double *input, double *output) {
  double bf0[32] = {0}, bf1[32] = {0};
  // stage 1;
  bf0[0] = input[0] / 4;
  bf0[1] = input[16] / 4;
  bf0[2] = input[8] / 4;
  bf0[3] = input[24] / 4;
  bf0[4] = input[4] / 4;
  bf0[5] = input[20] / 4;
  bf0[6] = input[12] / 4;
  bf0[7] = input[28] / 4;
  bf0[8] = input[2] / 4;
  bf0[9] = input[18] / 4;
  bf0[10] = input[10] / 4;
  bf0[11] = input[26] / 4;
  bf0[12] = input[6] / 4;
  bf0[13] = input[22] / 4;
  bf0[14] = input[14] / 4;
  bf0[15] = input[30] / 4;
  bf0[16] = input[1] / 4;
  bf0[17] = input[17] / 4;
  bf0[18] = input[9] / 4;
  bf0[19] = input[25] / 4;
  bf0[20] = input[5] / 4;
  bf0[21] = input[21] / 4;
  bf0[22] = input[13] / 4;
  bf0[23] = input[29] / 4;
  bf0[24] = input[3] / 4;
  bf0[25] = input[19] / 4;
  bf0[26] = input[11] / 4;
  bf0[27] = input[27] / 4;
  bf0[28] = input[7] / 4;
  bf0[29] = input[23] / 4;
  bf0[30] = input[15] / 4;
  bf0[31] = input[31] / 4;

  // stage 2
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = cospi[62] * bf0[16] - cospi[2] * bf0[31];
  bf1[17] = cospi[30] * bf0[17] - cospi[34] * bf0[30];
  bf1[18] = cospi[46] * bf0[18] - cospi[18] * bf0[29];
  bf1[19] = cospi[14] * bf0[19] - cospi[50] * bf0[28];
  bf1[20] = cospi[54] * bf0[20] - cospi[10] * bf0[27];
  bf1[21] = cospi[22] * bf0[21] - cospi[42] * bf0[26];
  bf1[22] = cospi[38] * bf0[22] - cospi[26] * bf0[25];
  bf1[23] = cospi[6] * bf0[23] - cospi[58] * bf0[24];
  bf1[24] = cospi[58] * bf0[23] + cospi[6] * bf0[24];
  bf1[25] = cospi[26] * bf0[22] + cospi[38] * bf0[25];
  bf1[26] = cospi[42] * bf0[21] + cospi[22] * bf0[26];
  bf1[27] = cospi[10] * bf0[20] + cospi[54] * bf0[27];
  bf1[28] = cospi[50] * bf0[19] + cospi[14] * bf0[28];
  bf1[29] = cospi[18] * bf0[18] + cospi[46] * bf0[29];
  bf1[30] = cospi[34] * bf0[17] + cospi[30] * bf0[30];
  bf1[31] = cospi[2] * bf0[16] + cospi[62] * bf0[31];

  // stage 3
  bf0[0] = bf1[0];
  bf0[1] = bf1[1];
  bf0[2] = bf1[2];
  bf0[3] = bf1[3];
  bf0[4] = bf1[4];
  bf0[5] = bf1[5];
  bf0[6] = bf1[6];
  bf0[7] = bf1[7];
  bf0[8] = cospi[60] * bf1[8] - cospi[4] * bf1[15];
  bf0[9] = cospi[28] * bf1[9] - cospi[36] * bf1[14];
  bf0[10] = cospi[44] * bf1[10] - cospi[20] * bf1[13];
  bf0[11] = cospi[12] * bf1[11] - cospi[52] * bf1[12];
  bf0[12] = cospi[52] * bf1[11] + cospi[12] * bf1[12];
  bf0[13] = cospi[20] * bf1[10] + cospi[44] * bf1[13];
  bf0[14] = cospi[36] * bf1[9] + cospi[28] * bf1[14];
  bf0[15] = cospi[4] * bf1[8] + cospi[60] * bf1[15];
  bf0[16] = bf1[16] + bf1[17];
  bf0[17] = bf1[16] - bf1[17];
  bf0[18] = -bf1[18] + bf1[19];
  bf0[19] = bf1[18] + bf1[19];
  bf0[20] = bf1[20] + bf1[21];
  bf0[21] = bf1[20] - bf1[21];
  bf0[22] = -bf1[22] + bf1[23];
  bf0[23] = bf1[22] + bf1[23];
  bf0[24] = bf1[24] + bf1[25];
  bf0[25] = bf1[24] - bf1[25];
  bf0[26] = -bf1[26] + bf1[27];
  bf0[27] = bf1[26] + bf1[27];
  bf0[28] = bf1[28] + bf1[29];
  bf0[29] = bf1[28] - bf1[29];
  bf0[30] = -bf1[30] + bf1[31];
  bf0[31] = bf1[30] + bf1[31];

  // stage 4
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = cospi[56] * bf0[4] - cospi[8] * bf0[7];
  bf1[5] = cospi[24] * bf0[5] - cospi[40] * bf0[6];
  bf1[6] = cospi[40] * bf0[5] + cospi[24] * bf0[6];
  bf1[7] = cospi[8] * bf0[4] + cospi[56] * bf0[7];
  bf1[8] = bf0[8] + bf0[9];
  bf1[9] = bf0[8] - bf0[9];
  bf1[10] = -bf0[10] + bf0[11];
  bf1[11] = bf0[10] + bf0[11];
  bf1[12] = bf0[12] + bf0[13];
  bf1[13] = bf0[12] - bf0[13];
  bf1[14] = -bf0[14] + bf0[15];
  bf1[15] = bf0[14] + bf0[15];
  bf1[16] = bf0[16];
  bf1[17] = -cospi[8] * bf0[17] + cospi[56] * bf0[30];
  bf1[18] = -cospi[56] * bf0[18] - cospi[8] * bf0[29];
  bf1[19] = bf0[19];
  bf1[20] = bf0[20];
  bf1[21] = -cospi[40] * bf0[21] + cospi[24] * bf0[26];
  bf1[22] = -cospi[24] * bf0[22] - cospi[40] * bf0[25];
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = -cospi[40] * bf0[22] + cospi[24] * bf0[25];
  bf1[26] = cospi[24] * bf0[21] + cospi[40] * bf0[26];
  bf1[27] = bf0[27];
  bf1[28] = bf0[28];
  bf1[29] = -cospi[8] * bf0[18] + cospi[56] * bf0[29];
  bf1[30] = cospi[56] * bf0[17] + cospi[8] * bf0[30];
  bf1[31] = bf0[31];

  // stage 5
  bf0[0] = cospi[32] * bf1[0] + cospi[32] * bf1[1];
  bf0[1] = cospi[32] * bf1[0] - cospi[32] * bf1[1];
  bf0[2] = cospi[48] * bf1[2] - cospi[16] * bf1[3];
  bf0[3] = cospi[16] * bf1[2] + cospi[48] * bf1[3];
  bf0[4] = bf1[4] + bf1[5];
  bf0[5] = bf1[4] - bf1[5];
  bf0[6] = -bf1[6] + bf1[7];
  bf0[7] = bf1[6] + bf1[7];
  bf0[8] = bf1[8];
  bf0[9] = -cospi[16] * bf1[9] + cospi[48] * bf1[14];
  bf0[10] = -cospi[48] * bf1[10] - cospi[16] * bf1[13];
  bf0[11] = bf1[11];
  bf0[12] = bf1[12];
  bf0[13] = -cospi[16] * bf1[10] + cospi[48] * bf1[13];
  bf0[14] = cospi[48] * bf1[9] + cospi[16] * bf1[14];
  bf0[15] = bf1[15];
  bf0[16] = bf1[16] + bf1[19];
  bf0[17] = bf1[17] + bf1[18];
  bf0[18] = bf1[17] - bf1[18];
  bf0[19] = bf1[16] - bf1[19];
  bf0[20] = -bf1[20] + bf1[23];
  bf0[21] = -bf1[21] + bf1[22];
  bf0[22] = bf1[21] + bf1[22];
  bf0[23] = bf1[20] + bf1[23];
  bf0[24] = bf1[24] + bf1[27];
  bf0[25] = bf1[25] + bf1[26];
  bf0[26] = bf1[25] - bf1[26];
  bf0[27] = bf1[24] - bf1[27];
  bf0[28] = -bf1[28] + bf1[31];
  bf0[29] = -bf1[29] + bf1[30];
  bf0[30] = bf1[29] + bf1[30];
  bf0[31] = bf1[28] + bf1[31];

  // stage 6
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = bf0[1] - bf0[2];
  bf1[3] = bf0[0] - bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = -cospi[32] * bf0[5] + cospi[32] * bf0[6];
  bf1[6] = cospi[32] * bf0[5] + cospi[32] * bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8] + bf0[11];
  bf1[9] = bf0[9] + bf0[10];
  bf1[10] = bf0[9] - bf0[10];
  bf1[11] = bf0[8] - bf0[11];
  bf1[12] = -bf0[12] + bf0[15];
  bf1[13] = -bf0[13] + bf0[14];
  bf1[14] = bf0[13] + bf0[14];
  bf1[15] = bf0[12] + bf0[15];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = -cospi[16] * bf0[18] + cospi[48] * bf0[29];
  bf1[19] = -cospi[16] * bf0[19] + cospi[48] * bf0[28];
  bf1[20] = -cospi[48] * bf0[20] - cospi[16] * bf0[27];
  bf1[21] = -cospi[48] * bf0[21] - cospi[16] * bf0[26];
  bf1[22] = bf0[22];
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = -cospi[16] * bf0[21] + cospi[48] * bf0[26];
  bf1[27] = -cospi[16] * bf0[20] + cospi[48] * bf0[27];
  bf1[28] = cospi[48] * bf0[19] + cospi[16] * bf0[28];
  bf1[29] = cospi[48] * bf0[18] + cospi[16] * bf0[29];
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 7
  bf0[0] = bf1[0] + bf1[7];
  bf0[1] = bf1[1] + bf1[6];
  bf0[2] = bf1[2] + bf1[5];
  bf0[3] = bf1[3] + bf1[4];
  bf0[4] = bf1[3] - bf1[4];
  bf0[5] = bf1[2] - bf1[5];
  bf0[6] = bf1[1] - bf1[6];
  bf0[7] = bf1[0] - bf1[7];
  bf0[8] = bf1[8];
  bf0[9] = bf1[9];
  bf0[10] = -cospi[32] * bf1[10] + cospi[32] * bf1[13];
  bf0[11] = -cospi[32] * bf1[11] + cospi[32] * bf1[12];
  bf0[12] = cospi[32] * bf1[11] + cospi[32] * bf1[12];
  bf0[13] = cospi[32] * bf1[10] + cospi[32] * bf1[13];
  bf0[14] = bf1[14];
  bf0[15] = bf1[15];
  bf0[16] = bf1[16] + bf1[23];
  bf0[17] = bf1[17] + bf1[22];
  bf0[18] = bf1[18] + bf1[21];
  bf0[19] = bf1[19] + bf1[20];
  bf0[20] = bf1[19] - bf1[20];
  bf0[21] = bf1[18] - bf1[21];
  bf0[22] = bf1[17] - bf1[22];
  bf0[23] = bf1[16] - bf1[23];
  bf0[24] = -bf1[24] + bf1[31];
  bf0[25] = -bf1[25] + bf1[30];
  bf0[26] = -bf1[26] + bf1[29];
  bf0[27] = -bf1[27] + bf1[28];
  bf0[28] = bf1[27] + bf1[28];
  bf0[29] = bf1[26] + bf1[29];
  bf0[30] = bf1[25] + bf1[30];
  bf0[31] = bf1[24] + bf1[31];

  // stage 8
  bf1[0] = bf0[0] + bf0[15];
  bf1[1] = bf0[1] + bf0[14];
  bf1[2] = bf0[2] + bf0[13];
  bf1[3] = bf0[3] + bf0[12];
  bf1[4] = bf0[4] + bf0[11];
  bf1[5] = bf0[5] + bf0[10];
  bf1[6] = bf0[6] + bf0[9];
  bf1[7] = bf0[7] + bf0[8];
  bf1[8] = bf0[7] - bf0[8];
  bf1[9] = bf0[6] - bf0[9];
  bf1[10] = bf0[5] - bf0[10];
  bf1[11] = bf0[4] - bf0[11];
  bf1[12] = bf0[3] - bf0[12];
  bf1[13] = bf0[2] - bf0[13];
  bf1[14] = bf0[1] - bf0[14];
  bf1[15] = bf0[0] - bf0[15];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = -cospi[32] * bf0[20] + cospi[32] * bf0[27];
  bf1[21] = -cospi[32] * bf0[21] + cospi[32] * bf0[26];
  bf1[22] = -cospi[32] * bf0[22] + cospi[32] * bf0[25];
  bf1[23] = -cospi[32] * bf0[23] + cospi[32] * bf0[24];
  bf1[24] = cospi[32] * bf0[23] + cospi[32] * bf0[24];
  bf1[25] = cospi[32] * bf0[22] + cospi[32] * bf0[25];
  bf1[26] = cospi[32] * bf0[21] + cospi[32] * bf0[26];
  bf1[27] = cospi[32] * bf0[20] + cospi[32] * bf0[27];
  bf1[28] = bf0[28];
  bf1[29] = bf0[29];
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 9
  output[0] = bf1[0] + bf1[31];
  output[1] = bf1[1] + bf1[30];
  output[2] = bf1[2] + bf1[29];
  output[3] = bf1[3] + bf1[28];
  output[4] = bf1[4] + bf1[27];
  output[5] = bf1[5] + bf1[26];
  output[6] = bf1[6] + bf1[25];
  output[7] = bf1[7] + bf1[24];
  output[8] = bf1[8] + bf1[23];
  output[9] = bf1[9] + bf1[22];
  output[10] = bf1[10] + bf1[21];
  output[11] = bf1[11] + bf1[20];
  output[12] = bf1[12] + bf1[19];
  output[13] = bf1[13] + bf1[18];
  output[14] = bf1[14] + bf1[17];
  output[15] = bf1[15] + bf1[16];
  output[16] = bf1[15] - bf1[16];
  output[17] = bf1[14] - bf1[17];
  output[18] = bf1[13] - bf1[18];
  output[19] = bf1[12] - bf1[19];
  output[20] = bf1[11] - bf1[20];
  output[21] = bf1[10] - bf1[21];
  output[22] = bf1[9] - bf1[22];
  output[23] = bf1[8] - bf1[23];
  output[24] = bf1[7] - bf1[24];
  output[25] = bf1[6] - bf1[25];
  output[26] = bf1[5] - bf1[26];
  output[27] = bf1[4] - bf1[27];
  output[28] = bf1[3] - bf1[28];
  output[29] = bf1[2] - bf1[29];
  output[30] = bf1[1] - bf1[30];
  output[31] = bf1[0] - bf1[31];
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

void exact_filter_32(const double *input, double *output, const double *h) {
  double temp[32] = {0};
  dct32(input, temp);
  for (int i = 0; i < 32; i++)
    temp[i] *= h[i];
  idct32(temp, output);
  return;
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
