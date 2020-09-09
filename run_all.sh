#!/bin/bash

echo "Tikhonov, 4x4 grid..."
./speed_tikhonov4x4 input/data16_20000.txt output/rt_tik4x4.txt
echo "Tikhonov, 8x8 grid..."
./speed_tikhonov8x8 input/data64_20000.txt output/rt_tik8x8.txt
echo "Tikhonov, 16x16 grid..."
./speed_tikhonov16x16 input/data256_20000.txt output/rt_tik16x16.txt
echo "Tikhonov, length-32 line..."
./speed_tikhonov32 input/data32_20000.txt output/rt_tik32.txt
echo "Tikhonov, length-64 line..."
./speed_tikhonov64 input/data64_20000.txt output/rt_tik64.txt
echo "Tikhonov, length-128 line..."
./speed_tikhonov128 input/data128_20000.txt output/rt_tik128.txt

echo "Diffusion, 4x4 grid..."
./speed_diffusion4x4 input/data16_20000.txt output/rt_diff4x4.txt
echo "Diffusion, 8x8 grid..."
./speed_diffusion8x8 input/data64_20000.txt output/rt_diff8x8.txt
echo "Diffusion, 16x16 grid..."
./speed_diffusion16x16 input/data64_20000.txt output/rt_diff16x16.txt
echo "Diffusion, length-32 line..."
./speed_diffusion32 input/data32_20000.txt output/rt_diff32.txt
echo "Diffusion, length-64 line..."
./speed_diffusion64 input/data64_20000.txt output/rt_diff64.txt
echo "Diffusion, length-128 line..."
./speed_diffusion128 input/data128_20000.txt output/rt_diff128.txt

echo "Ideal low-pass, 4x4 grid..."
./speed_ideallowpass4x4 input/data16_20000.txt output/rt_lp4x4.txt
echo "Ideal low-pass, 8x8 grid..."
./speed_ideallowpass8x8 input/data64_20000.txt output/rt_lp8x8.txt
echo "Ideal low-pass, 16x16 grid..."
./speed_ideallowpass16x16 input/data64_20000.txt output/rt_lp16x16.txt
echo "Ideal low-pass, length-32 line..."
./speed_ideallowpass32 input/data32_20000.txt output/rt_lp32.txt
echo "Ideal low-pass, length-64 line..."
./speed_ideallowpass64 input/data64_20000.txt output/rt_lp64.txt
echo "Ideal low-pass, length-128 line..."
./speed_ideallowpass128 input/data128_20000.txt output/rt_lp128.txt

echo "Exponential, 4x4 grid..."
./speed_exp4x4 input/data16_20000.txt output/rt_exp4x4.txt
echo "Exponential, 8x8 grid..."
./speed_exp8x8 input/data64_20000.txt output/rt_exp8x8.txt
echo "Exponential, 16x16 grid..."
./speed_exp16x16 input/data64_20000.txt output/rt_exp16x16.txt
echo "Exponential, length-32 line..."
./speed_exp32 input/data32_20000.txt output/rt_exp32.txt
echo "Exponential, length-64 line..."
./speed_exp64 input/data64_20000.txt output/rt_exp64.txt
echo "Exponential, length-128 line..."
./speed_exp128 input/data128_20000.txt output/rt_exp128.txt
