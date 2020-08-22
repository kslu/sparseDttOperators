#include "../src/grfilter.h"

#define LEN 16

int main(int argc, char *argv[]) {

  double input4[4] = {0.7536, 0.8593, 0.7907, 0.9612};
  double input8[8] = {0.2638, 0.5868, 0.5683, 0.7702,
                      0.3364, 0.7005, 0.1097, 0.2274};
  double input16[16] = {0.8537, 0.4644, 0.5917, 0.3501, 0.1662, 0.6533,
                        0.9051, 0.1245, 0.2605, 0.9276, 0.8231, 0.7579,
                        0.6901, 0.3326, 0.8131, 0.3530};
  double input32[32] = {0.3041, 0.4868, 0.9101, 0.8521, 0.9781, 0.6740, 0.9982,
                        0.1307, 0.7803, 0.9581, 0.6719, 0.4030, 0.6437, 0.0349,
                        0.7468, 0.2155, 0.8643, 0.6717, 0.6606, 0.9273, 0.1661,
                        0.1734, 0.6295, 0.3849, 0.9915, 0.6498, 0.5767, 0.3216,
                        0.8256, 0.0776, 0.2773, 0.6198};
  double input64[64] = {
      0.1523, 0.0659, 0.8370, 0.6498, 0.8000, 0.0484, 0.9789, 0.5670,
      0.0113, 0.2244, 0.6273, 0.7424, 0.1177, 0.0295, 0.8463, 0.0485,
      0.6853, 0.3698, 0.9307, 0.5531, 0.9579, 0.5442, 0.8037, 0.2446,
      0.4558, 0.1768, 0.4685, 0.4764, 0.5873, 0.3288, 0.9122, 0.5572,
      0.3888, 0.7945, 0.1142, 0.7704, 0.5856, 0.5918, 0.5072, 0.1513,
      0.5907, 0.9409, 0.7874, 0.5198, 0.4786, 0.5121, 0.2807, 0.8433,
      0.1618, 0.8253, 0.4243, 0.6659, 0.1212, 0.9474, 0.4233, 0.6638,
      0.6048, 0.1414, 0.3167, 0.4626, 0.3348, 0.4912, 0.9928, 0.5254};
  double output4[4] = {0}, output8[8] = {0}, output16[16] = {0},
         output32[32] = {0}, output64[64] = {0};
  double recon4[4] = {0}, recon8[8] = {0}, recon16[16] = {0}, recon32[32] = {0},
         recon64[64] = {0};

  // test dct4
  fprintf(stderr, "=== dct4 ===\nInput: [");
  for (int i = 0; i < 4; i++)
    fprintf(stderr, "%f, ", input4[i]);
  fprintf(stderr, "]\n");
  dct4(input4, output4);
  // MATLAB output: [ 1.6824, -0.1170, 0.0324, -0.1010]
  fprintf(stderr, "Output: [");
  for (int i = 0; i < 4; i++)
    fprintf(stderr, "%f, ", output4[i]);
  fprintf(stderr, "]\n");
  idct4(output4, recon4);
  fprintf(stderr, "Recon: [");
  for (int i = 0; i < 4; i++)
    fprintf(stderr, "%f, ", recon4[i]);
  fprintf(stderr, "]\n");

  // test dct8
  fprintf(stderr, "=== dct8 ===\nInput: [");
  for (int i = 0; i < 8; i++)
    fprintf(stderr, "%f, ", input8[i]);
  fprintf(stderr, "]\n");
  dct8(input8, output8);
  // MATLAB output: [ 1.2597, 0.2218, -0.3938, -0.0871, -0.1300, -0.0564,
  //     0.1466, -0.3967];
  fprintf(stderr, "Output: [");
  for (int i = 0; i < 8; i++)
    fprintf(stderr, "%f, ", output8[i]);
  fprintf(stderr, "]\n");
  idct8(output8, recon8);
  fprintf(stderr, "Recon: [");
  for (int i = 0; i < 8; i++)
    fprintf(stderr, "%f, ", recon8[i]);
  fprintf(stderr, "]\n");

  // test dct4x4
  fprintf(stderr, "=== dct4x4 ===\nInput: [");
  for (int i = 0; i < 16; i++)
    fprintf(stderr, "%f, ", input16[i]);
  fprintf(stderr, "]\n");
  dct4x4(input16, output16);
  // MATLAB output: [ 2.2668, 0.0236, -0.4887, 0.2988, -0.1013, 0.2027, 0.0093,
  //     -0.0188, -0.0424, 0.3612, 0.5113, 0.3259, 0.3101, -0.1118, 0.2087,
  //     -0.2976 ]
  fprintf(stderr, "Output: [");
  for (int i = 0; i < 16; i++)
    fprintf(stderr, "%f, ", output16[i]);
  fprintf(stderr, "]\n");
  idct4x4(output16, recon16);
  fprintf(stderr, "Recon: [");
  for (int i = 0; i < 16; i++)
    fprintf(stderr, "%f, ", recon16[i]);
  fprintf(stderr, "]\n");

  // test dct8x8
  fprintf(stderr, "=== dct8x8 ===\nInput: [");
  for (int i = 0; i < 64; i++)
    fprintf(stderr, "%f, ", input64[i]);
  fprintf(stderr, "]\n");
  dct8x8(input64, output64);
  // MATLAB output: [ 4.0949, -0.2907, -0.2677, -0.2645, -0.2260, 0.4894,
  // -0.3567, 0.1979, -0.1779, -0.0422, -0.4884, -0.3040, 0.0721, 0.4463,
  // -0.2292, 0.8029, -0.1972, -0.3482, 0.0346, -0.3105, -0.0790, 0.3763,
  // 0.0454, 0.1719, 0.0591, 0.1639, -0.0204, -0.3749, 0.0213, -0.4967, -0.0476,
  // -0.3570, -0.1359, -0.3710, 0.0406, 0.2528, 0.3018, 0.2728, -0.2629, 0.2862,
  // 0.3109, -0.2662, 0.0671, -0.0919, 0.0305, -0.2178, 0.2756, 0.1029, 0.5234,
  // 0.1566, 0.0990, 0.0862, 0.2579, -0.1671, 0.2502, 0.5815, 0.1733, 0.1306,
  // -0.3298, 0.3408, -0.0807, -0.1660, -0.0775, -0.2163 ]
  fprintf(stderr, "Output: [");
  for (int i = 0; i < 64; i++)
    fprintf(stderr, "%f, ", output64[i]);
  fprintf(stderr, "]\n");
  idct8x8(output64, recon64);
  fprintf(stderr, "Recon: [");
  for (int i = 0; i < 64; i++)
    fprintf(stderr, "%f, ", recon64[i]);
  fprintf(stderr, "]\n");

  // test dct32
  fprintf(stderr, "=== dct32 ===\nInput: [");
  for (int i = 0; i < 32; i++)
    fprintf(stderr, "%f, ", input32[i]);
  fprintf(stderr, "]\n");
  dct32(input32, output32);
  // MATLAB output: [ 3.2891, 0.3984, 0.0567, 0.1074, -0.2785, -0.1890, -0.4170,
  // 0.1018, -0.0341, -0.7801, -0.1111, -0.0611, -0.0980, 0.1544, 0.1844,
  // -0.2942, 0.0374, -0.3367, 0.0088, -0.1898, 0.4855, 0.0357, 0.1131, 0.2220,
  // -0.4577, -0.3398, 0.2405, 0.2958, -0.2917, 0.2847, 0.1382, 0.6170 ]
  fprintf(stderr, "Output: [");
  for (int i = 0; i < 32; i++)
    fprintf(stderr, "%f, ", output32[i]);
  fprintf(stderr, "]\n");
  idct32(output32, recon32);
  fprintf(stderr, "Recon: [");
  for (int i = 0; i < 32; i++)
    fprintf(stderr, "%f, ", recon32[i]);
  fprintf(stderr, "]\n");

  return 0;
}
