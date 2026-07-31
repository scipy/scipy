#!/bin/bash
# Run the SciPy test suite for the WebAssembly/Pyodide build, in splits to
# avoid memory corruption issues noticed when not doing so, and combine the
# logs, see gh-25713

logdir="${RUNNER_TEMP:-${TMPDIR:-/tmp}}/scipy-reportlogs"
mkdir -p "$logdir"
ret=0

echo "Writing per-group report logs to $logdir"

# number of pytest-split groups to run
N=10

# iteration number of the current group
i=1

while [ "$i" -le "$N" ]; do
  echo "::group::SciPy test group $i of $N"
  python -m pytest -vra --pyargs scipy -m 'not slow' --splits "$N" --group "$i" --report-log="$logdir/group-$i.jsonl"
  rc=$?
  echo "::endgroup::"
  if [ "$rc" -ne 0 ] && [ "$rc" -ne 5 ]; then ret=1; echo "!! pytest-split group $i failed (exit $rc)"; fi
  i=$((i + 1))
done

echo "Combined SciPy Pyodide test suite report log: $logdir/full.jsonl"
cat "$logdir"/group-*.jsonl > "$logdir/full.jsonl" 2>/dev/null || true

exit $ret
