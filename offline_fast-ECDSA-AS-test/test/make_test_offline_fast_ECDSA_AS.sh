#!/bin/bash
dir="~/Documents/Project/openssl-master"

g++ -std=c++11 -pthread -O3 test_offline_fast_ECDSA_AS.cpp -L ${dir} -l ssl -l crypto -o test_fast_offline_ECDSA_AS -I ${dir}