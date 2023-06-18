set -x

bsub -b -I -q q_sw_cpc2023 -shared -n 1 -cgsp 64 -host_stack 1024 -share_size 15000 -cache_size 32 ./pcg_solve
