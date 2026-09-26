# Benchmarks

`discrete_log` first computes the order of the base (which factors the modulus and every
`p - 1`), then picks an algorithm from the order. Either part can dominate, so they are
measured separately:

| What | Where | Measure | CI |
|---|---|---|---|
| Building blocks: factorization (trial division, Pollard's rho, perfect powers), `n_order`, order of the base in one pass, index calculus smoothness test | `ops.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Each algorithm alone, the order given: the three prime-order algorithms on the same instances from 28 to 46 bits, index calculus, Pohlig-Hellman | `ops.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Complete `discrete_log` on the instances of the tests | `e2e.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Complete `discrete_log` on larger instances | `walltime.rs` | Wall-clock time (criterion) | ❌ too noisy on shared runners |

Pollard's rho and index calculus are randomized: the benchmarks use their seeded variants, so
every run draws the same numbers, and their instances solve several targets to average out the
luck of one draw.

The three prime-order algorithms run on the same safe-prime instances at 28, 32, 34, 36, 40, 42
and 46 bits, which is how the selection boundaries in `discrete_log` were chosen (baby-step
giant-step stops at 40 bits: above that the order passes its own memory cap and it refuses the
problem). `element_order power_of_two` covers the exact Carmichael lambda of a modulus divisible
by 8, and `pohlig_hellman prime_power_order` (base 2 modulo 1009^3) covers the baby-step table
shared between the digits of a repeated prime.

`[profile.bench]` sets `codegen-units = 1`: without it, an edit in one module shifts the one-time
costs of another by up to 2% and the comparison reports regressions on benchmarks that ran no new
code. It costs a few seconds per benchmark rebuild.

## Instruction counts

Requires [Valgrind](https://valgrind.org/) and `gungraun-runner` with the same version as the
`gungraun` dev-dependency:

```sh
cargo install gungraun-runner --version 0.19.4 --locked

cargo bench --features bench --bench ops --bench e2e -- --parallel=auto
# Only some benchmarks: FILE::GROUP::FUNCTION::ID wildcard
cargo bench --features bench --bench ops -- 'ops::prime_order::*'
```

Compare two versions of the code:

```sh
git checkout main && cargo bench --features bench --bench ops --bench e2e -- --save-baseline=main
git checkout -    && cargo bench --features bench --bench ops --bench e2e -- --baseline=main
```

## Wall-clock time

Instruction counts don't see cache effects nor the build profile: inlining across crates with
`lto = "fat"`, less register pressure with `codegen-units = 1`. `walltime.rs` measures the
largest instances, with the default profile or with the `bench-lto` profile (`bench` +
`lto = "fat"`, `codegen-units = 1`):

```sh
cargo bench --features bench --bench walltime -- --save-baseline before
cargo bench --features bench --bench walltime -- --baseline before
cargo bench --profile bench-lto --features bench --bench walltime
```

## CI

`.github/workflows/bench.yml`:

- On pull requests, the base branch and the PR run in the same job and are compared there. The
  job fails if an instruction count grows by more than 2%.
- On pull requests, one comment compares the PR with its base branch, worst changes first and
  unchanged benchmarks collapsed (updated on each push).
- Results are published whatever happens. A regression, a benchmark that fails outright and a
  baseline that could not be measured all leave the numbers in the job summary, in the comment and
  on the history charts, and the job still goes red afterwards. A failed run says so in the
  comment, since its numbers are only what it managed to measure.
- The benchmarks are built without debug info in CI (`CARGO_PROFILE_BENCH_DEBUG=false`): half
  the build time, same instruction counts.
- On pushes to `main`, results are stored on the `gh-pages` branch by
  [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark).
