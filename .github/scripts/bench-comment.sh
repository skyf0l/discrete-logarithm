#!/usr/bin/env bash
# Builds the benchmark comparison comment for the PR (markdown on stdout), from the downloaded
# artifact instructions/instructions.md (gungraun-comment.sh).
#
# Env: RUN_URL (link to the workflow run), BENCH_FAILED (a benchmark failed, so the numbers below
# are what the run managed to measure).
set -euo pipefail

echo '## Benchmarks: PR vs base branch'
echo
if [ "${BENCH_FAILED:-false}" = true ]; then
  echo '> [!WARNING]'
  echo '> Some benchmarks failed. What follows is what the run measured, which may be partial.'
  echo
fi
cat instructions/instructions.md 2>/dev/null || echo '_No instruction count comparison._'
echo
echo "<sub>[Full results]($RUN_URL) · Instruction counts are exact (Valgrind); estimated" \
  "cycles weight memory accesses with a simulated cache.</sub>"
