//! Wall-clock benchmarks (criterion) of complete discrete logarithms, for local use.
//!
//! Complements the instruction-count benchmarks with effects Valgrind does not measure
//! (cache, allocator, CPU-specific GMP code), on instances too large to run under Valgrind.
//! Too noisy for shared CI runners: compare locally with
//! `cargo bench --features bench --bench walltime -- --save-baseline before` then `-- --baseline before`.

use criterion::{criterion_group, criterion_main, Criterion};
use discrete_logarithm::discrete_log;
use rug::Integer;
use std::{hint::black_box, str::FromStr};

/// `(name, n, a, b)`, solved with `discrete_log`.
const INSTANCES: [(&str, &str, &str, &str); 2] = [
    // Large prime order (about 2^47): Pollard's rho.
    (
        "large_prime_order",
        "265390227570863",
        "184500076053622",
        "2",
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
