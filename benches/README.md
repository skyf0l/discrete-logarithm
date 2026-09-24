# Benchmarks

`discrete_log` first computes the order of the base (which factors the modulus and every
`p - 1`), then picks an algorithm from the order. Either part can dominate, so they are
measured separately:

| What | Where | Measure | CI |
|---|---|---|---|
| Building blocks: factorization, `n_order`, index calculus smoothness test | `ops.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Each algorithm alone, the order given: the prime-order algorithms on the same instances, index calculus, Pohlig-Hellman | `ops.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Complete `discrete_log` on the instances of the tests | `e2e.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Complete `discrete_log` on larger instances | `walltime.rs` | Wall-clock time (criterion) | ❌ too noisy on shared runners |

Pollard's rho and index calculus are randomized but always draw the same random numbers, so
their instances solve several targets to average out the luck.

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
`lto = "fat"`, less register pressure with `codegen-units = 1`. `walltime.rs` measures the two
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
  unchanged benchmarks collapsed (updated on each push). A regression fails the job but still
  reports the results.
- The benchmarks are built without debug info in CI (`CARGO_PROFILE_BENCH_DEBUG=false`): half
  the build time, same instruction counts.
- On pushes to `main`, results are stored on the `gh-pages` branch by
  [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark).
