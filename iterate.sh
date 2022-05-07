#!/bin/bash
if [[ "$#" -lt 3 ]]; then
  echo "Usage: iterate <start> <diff> <count> [rest of params to be passed to main]"
  echo "--tempAmb <actual value> will be appended."
  echo "Do not specify output name, it will be series<n>.png"
  exit
fi
i=0
t=$1
d=$2
n=$3
shift
shift
shift
while [[ $i -lt $n ]]; do
  ./main --tempAmb $t $*
  mv result.png series$(printf "%03d" $i).png
  i=$((i+1))
  t=`awk "BEGIN{print $t + $d}"`
done
