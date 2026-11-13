#!/bin/bash
# Run "target title <job index>" njob times, ncore at a time (simple version of runmulti.py).
# Usage: ./runmulti.sh target title njob [ncore]
if [ $# -lt 3 ]; then
  echo "Usage: $0 target title njob [ncore]" >&2
  exit 2
fi
target=$1; title=$2; njob=$3; ncore=${4:-4}
echo RUN $target $title : $njob jobs, $ncore at a time

trap 'kill $(jobs -p) 2>/dev/null; exit 130' INT TERM

fail=0
pids=""
for (( i = 0; i < njob; i++ )); do
  while [ "$(jobs -rp | wc -l)" -ge "$ncore" ]; do
    sleep 0.5
  done
  $target $title $i &
  pids="$pids $!"
done
# Wait by the recorded pids: "jobs -p" no longer lists the jobs that already
# finished, so their exit codes would be lost.
for p in $pids; do
  wait $p || fail=1
done
exit $fail
