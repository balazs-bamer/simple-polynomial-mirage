#!/bin/bash
if [[ "$#" -lt 4 ]]; then
  echo "Usage: bash iterateMain.sh <start> <diff> <count> <parameterToIterate> [rest of params to be passed to main]"
  echo "parameterToIterate <actual value> will be appended."
  echo "Do not specify output name, it will be series<n>.png"
  exit
fi
i=0
t=$1
d=$2
n=$3
p=$4
shift
shift
shift
shift
while [[ $i -lt $n ]]; do
  i=$((i+1))
  echo Iteration $i of $n
  ./main $p $t $*
  mv result.png series$(printf "%03d" $i).png
  t=`awk "BEGIN{print $t + $d}" | tr ',' '.'`
done
