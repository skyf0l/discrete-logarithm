#!/usr/bin/env bash
# Builds the benchmark comparison comment for the PR (markdown on stdout), from the downloaded
# artifact instructions/instructions.md (gungraun-comment.sh).
#
# Env: RUN_URL (link to the workflow run).
set -euo pipefail

echo '## Benchmarks: PR vs base branch'
echo
cat instructions/instructions.md 2>/dev/null || echo '_No instruction count comparison._'
echo
echo "<sub>[Full results]($RUN_URL) · Instruction counts are exact (Valgrind); estimated" \
  "cycles weight memory accesses with a simulated cache.</sub>"
