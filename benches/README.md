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

```sh
cargo bench --features bench --bench walltime -- --save-baseline before
cargo bench --features bench --bench walltime -- --baseline before
```

## CI

`.github/workflows/bench.yml`:

- On pull requests, the base branch and the PR run in the same job and are compared there. The
  job fails if an instruction count grows by more than 2%.
- On pull requests, github-action-benchmark also comments a comparison with the last `main` results
  on the PR (updated on each push).
- On pushes to `main`, results are stored on the `gh-pages` branch by
  [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark).
