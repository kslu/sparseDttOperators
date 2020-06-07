#!/bin/bash

./speed_tikhonov4x4 input/data16_20000.txt output/rt_tik4x4.txt
./speed_tikhonov8x8 input/data64_20000.txt output/rt_tik8x8.txt
./speed_tikhonov32 input/data32_20000.txt output/rt_tik32.txt
./speed_tikhonov128 input/data128_20000.txt output/rt_tik128.txt
./speed_ideallowpass4x4 input/data16_20000.txt output/rt_lp4x4.txt
./speed_ideallowpass8x8 input/data64_20000.txt output/rt_lp8x8.txt
./speed_ideallowpass32 input/data32_20000.txt output/rt_lp32.txt
./speed_ideallowpass128 input/data128_20000.txt output/rt_lp128.txt
