//! Instruction-count benchmarks (Valgrind/Callgrind) of complete discrete logarithms.
//!
//! `discrete_log` computes the order of the base (factorization of the modulus and of the
//! `p - 1`), then chooses the algorithm: these benchmarks measure everything together, on the
//! instances of the tests. Results are checked: a wrong logarithm makes the benchmark fail
//! instead of reporting a fake speedup.
//!
//! Instances are small on purpose: Valgrind slows execution down by about 50x.
//!
//! Run with `cargo bench --features bench --bench e2e` (requires Valgrind and `gungraun-runner`).

mod common;

use common::{DIGITS_108, int};
use discrete_logarithm::discrete_log;
use gungraun::{library_benchmark, library_benchmark_group, main};
use rug::Integer;
use std::hint::black_box;

/// `(n, a, b)`: base `b` and target `b^x mod n`.
fn power(n: &str, b: u32, x: u32) -> (Integer, Integer, Integer) {
    let n = int(n);
    let b = Integer::from(b);
    let a = b.clone().pow_mod(&Integer::from(x), &n).unwrap();
    (n, a, b)
}

#[library_benchmark]
// Order below 1000: trial multiplication.
#[bench::n_587(power("587", 2, 9))]
// Prime order: baby-step giant-step.
#[bench::n_2456747(power("2456747", 3, 51))]
// Composite modulus, composite order: Pohlig-Hellman.
#[bench::n_32942478(power("32942478", 11, 127))]
#[bench::n_5779((int("5779"), int("3528"), int("6215")))]
// Two prime factors above the trial division bound: Pollard's rho factors the modulus.
#[bench::two_large_primes((int("1125902456980891"), int("55959216458270"), int("3")))]
// Large modulus, smooth order: Pohlig-Hellman on large numbers.
#[bench::digits_108((int(DIGITS_108.0), int(DIGITS_108.1), int(DIGITS_108.2)))]
fn solve(input: (Integer, Integer, Integer)) -> Integer {
    let (n, a, b) = &input;
    let x = black_box(discrete_log(black_box(n), black_box(a), black_box(b))).unwrap();
    assert_eq!(
        &b.clone().pow_mod(&x, n).unwrap(),
        a,
        "wrong logarithm {x} of {a} in base {b} modulo {n}"
    );
    x
}

library_benchmark_group!(name = e2e, benchmarks = [solve]);

main!(library_benchmark_groups = e2e);
