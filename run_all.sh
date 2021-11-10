#!/bin/bash

echo "Tikhonov, 16x16 grid..."
./speed_tikhonov16x16 input/data256_20000.txt output/rt_tik16x16.txt
echo "Tikhonov, length-64 line..."
./speed_tikhonov64 input/data64_20000.txt output/rt_tik64.txt

echo "Ideal low-pass, 16x16 grid..."
./speed_ideallowpass16x16 input/data64_20000.txt output/rt_lp16x16.txt
echo "Ideal low-pass, length-64 line..."
./speed_ideallowpass64 input/data64_20000.txt output/rt_lp64.txt

echo "Exponential, 16x16 grid..."
./speed_exp16x16 input/data256_20000.txt output/rt_exp16x16.txt
echo "Exponential, length-64 line..."
./speed_exp64 input/data64_20000.txt output/rt_exp64.txt
