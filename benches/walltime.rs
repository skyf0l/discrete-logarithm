//! Wall-clock benchmarks (criterion) of complete discrete logarithms, for local use.
//!
//! Complements the instruction-count benchmarks with what Valgrind does not see: cache and
//! memory effects, CPU-specific GMP code, and the build profile (inlining across crates with
//! `lto = "fat"`), on instances too large to run under Valgrind. Too noisy for shared CI
//! runners.
//!
//! ```text
//! # Default bench profile
//! cargo bench --features bench --bench walltime -- --save-baseline before
//! cargo bench --features bench --bench walltime -- --baseline before
//! # With the optimizations of a final build (`lto = "fat"`, `codegen-units = 1`)
//! cargo bench --profile bench-lto --features bench --bench walltime
//! ```

use criterion::{Criterion, criterion_group, criterion_main};
use discrete_logarithm::discrete_log;
use rug::Integer;
use std::{hint::black_box, str::FromStr};

/// `(name, n, a, b)`, solved with `discrete_log`.
const INSTANCES: [(&str, &str, &str, &str); 3] = [
    // Large prime order (about 2^47): Pollard's rho.
    (
        "large_prime_order",
        "265390227570863",
        "184500076053622",
        "2",
    ),
    // Two prime factors above the trial division bound: Pollard's rho factors the modulus.
    (
        "two_large_primes",
        "1125902456980891",
        "55959216458270",
        "3",
    ),
    // Large modulus, smooth order: Pohlig-Hellman on large numbers.
    (
        "digits_108",
        "22708823198678103974314518195029102158525052496759285596453269189798311427475159776411276642277139650833937",
        "17463946429475485293747680247507700244427944625055089103624311227422110546803452417458985046168310373075327",
        "123456",
    ),
];

fn solve(c: &mut Criterion) {
    let mut group = c.benchmark_group("discrete_log");
    group.sample_size(10);

    for (name, n, a, b) in INSTANCES {
        let n = Integer::from_str(n).unwrap();
        let a = Integer::from_str(a).unwrap();
        let b = Integer::from_str(b).unwrap();
        group.bench_function(name, |bencher| {
            bencher.iter(|| {
                black_box(discrete_log(black_box(&n), black_box(&a), black_box(&b))).unwrap()
            })
        });
    }

    group.finish();
}

criterion_group!(benches, solve);
criterion_main!(benches);
