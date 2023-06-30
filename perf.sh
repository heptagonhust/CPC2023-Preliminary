#!/bin/bash

prefix="slave"
output_dir="perf.d"
output_file="perf.out"

cd $(dirname "$0") || exit 1

make
./run.sh | grep "^$prefix" | tr -d '\r' > "$output_file"
if [[ ! -d "$output_dir" ]]; then
  mkdir -p "$output_dir"
fi
for i in $(seq 0 63); do
  (
    grep "$(printf "slave%02d" "$i")" "$output_file" > "$(printf "%s/slave%02d" "$output_dir" "$i")"
  )
done
wait
