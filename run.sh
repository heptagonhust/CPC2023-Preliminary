set -x

bsub -b -I -q q_sw_cpc2023 -shared -n 1 -cgsp 64 -host_stack 1024 -share_size 15000 -ldm_share_mode 5 -ldm_share_size 16 ./pcg_solve
