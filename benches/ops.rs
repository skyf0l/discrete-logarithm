//! Instruction-count benchmarks (Valgrind/Callgrind) of the building blocks of `discrete_log`.
//!
//! Every input is deterministic (fixed seeds), and the setup work (generating instances,
//! computing the order given to the algorithms) is not measured. Results are checked: a wrong
//! answer makes the benchmark fail instead of reporting a fake speedup.
//!
//! Run with `cargo bench --features bench --bench ops` (requires Valgrind and `gungraun-runner`).

mod common;

use common::{check_factors, int, prime_bits, Instance, DIGITS_108, SEED};
use discrete_logarithm::{
    bench::{fast_factor, is_smooth},
    discrete_log_index_calculus, discrete_log_pohlig_hellman, discrete_log_pollard_rho,
    discrete_log_shanks_steps, n_order, Error,
};
use gungraun::{library_benchmark, library_benchmark_group, main};
use primal::Primes;
use rug::Integer;
use std::{collections::HashMap, hint::black_box};

/// Product of the first primes with small exponents: factored by trial division alone.
fn smooth() -> Integer {
    [
        (2, 10),
        (3, 5),
        (5, 3),
        (7, 2),
        (11, 1),
        (13, 1),
        (17, 1),
        (19, 1),
    ]
    .iter()
    .fold(Integer::from(1), |acc, &(p, e)| {
        acc * Integer::from(Integer::u_pow_u(p, e))
    })
}

/// Small factor times a prime larger than every trial divisor: the worst case of trial division.
fn large_prime_cofactor() -> Integer {
    prime_bits(40, SEED) * 587u32
}

#[library_benchmark]
#[bench::small(int("587"))]
#[bench::smooth(smooth())]
#[bench::large_prime_cofactor(large_prime_cofactor())]
fn factor(n: Integer) -> HashMap<Integer, usize> {
    let factors = black_box(fast_factor(black_box(&n)));
    check_factors(&n, &factors);
    factors
}

/// Factor base of the index calculus for `n = 999231337607`: the primes below 508.
fn factorbase() -> Vec<usize> {
    Primes::all().take_while(|&p| p < 508).collect()
}

/// A number of 39 bits (the size of the numbers tested), smooth over the factor base.
fn smooth_number() -> Integer {
    Integer::from(2u64.pow(3) * 3u64.pow(2) * 7 * 97 * 101 * 211 * 499)
}

/// A number of about 39 bits with a large prime factor: the common case, where every prime of the
/// factor base is tried.
fn not_smooth_number() -> Integer {
    prime_bits(36, SEED) * 6u32
}

#[library_benchmark]
#[bench::smooth((smooth_number(), factorbase()))]
#[bench::not_smooth((not_smooth_number(), factorbase()))]
fn index_calculus_smoothness(input: (Integer, Vec<usize>)) -> Option<Vec<u32>> {
    let (n, factorbase) = black_box(input);
    black_box(is_smooth(n, &factorbase))
}

library_benchmark_group!(
    name = factorization,
    benchmarks = [factor, index_calculus_smoothness]
);

#[library_benchmark]
#[bench::prime((int("3"), int("2456747")))]
#[bench::prime_power((int("2"), Integer::from(Integer::u_pow_u(3, 40))))]
#[bench::composite((int("11"), int("32942478")))]
fn order(input: (Integer, Integer)) -> Integer {
    let (a, n) = black_box(&input);
    black_box(n_order(a, n).unwrap())
}

library_benchmark_group!(name = n_order_group, benchmarks = [order]);

type Algorithm = fn(&Integer, &Integer, &Integer, Option<&Integer>) -> Result<Integer, Error>;

/// Solves every target of `instance` with `algorithm`, the order given, and checks the results.
fn solve(algorithm: Algorithm, instance: &Instance) -> Vec<Integer> {
    let logs: Vec<Integer> = instance
        .targets
        .iter()
        .map(|a| {
            black_box(algorithm(
                black_box(&instance.n),
                black_box(a),
                black_box(&instance.b),
                Some(&instance.order),
            ))
            .unwrap()
        })
        .collect();
    instance.check(&logs);
    logs
}

// The three algorithms for prime orders, on the same instances: shows where each one is the
// fastest, to choose the thresholds of `discrete_log_with_order`.

#[library_benchmark]
#[bench::bits_28(Instance::safe_prime(28, SEED))]
#[bench::bits_34(Instance::safe_prime(34, SEED))]
fn prime_order_shanks_steps(instance: Instance) -> Vec<Integer> {
    solve(discrete_log_shanks_steps, &instance)
}

#[library_benchmark]
#[bench::bits_28(Instance::safe_prime(28, SEED))]
#[bench::bits_34(Instance::safe_prime(34, SEED))]
fn prime_order_pollard_rho(instance: Instance) -> Vec<Integer> {
    solve(discrete_log_pollard_rho, &instance)
}

#[library_benchmark]
#[bench::bits_28(Instance::safe_prime(28, SEED))]
#[bench::bits_34(Instance::safe_prime(34, SEED))]
fn prime_order_index_calculus(instance: Instance) -> Vec<Integer> {
    solve(discrete_log_index_calculus, &instance)
}

// The index calculus instances of the tests.
#[library_benchmark]
#[bench::n_983(Instance::with_order("983", "948", "2", "491"))]
#[bench::n_633383(Instance::with_order("633383", "21794", "2", "316691"))]
#[bench::n_941762639(Instance::with_order("941762639", "68822582", "2", "470881319"))]
#[bench::n_999231337607(Instance::with_order("999231337607", "888188918786", "2", "499615668803"))]
#[bench::n_47747730623(Instance::with_order(
    "47747730623",
    "19410045286",
    "43425105668",
    "645239603"
))]
fn index_calculus(instance: Instance) -> Vec<Integer> {
    solve(discrete_log_index_calculus, &instance)
}

library_benchmark_group!(
    name = prime_order,
    benchmarks = [
        prime_order_shanks_steps,
        prime_order_pollard_rho,
        prime_order_index_calculus,
        index_calculus
    ]
);

// Composite orders, the order given: factorization of the order and one discrete logarithm
// per prime power.
#[library_benchmark]
#[bench::n_32942478(Instance::new("32942478", "11", &["10792037"]))]
#[bench::digits_108(Instance::new(DIGITS_108.0, DIGITS_108.2, &[DIGITS_108.1]))]
fn pohlig_hellman(instance: Instance) -> Vec<Integer> {
    solve(discrete_log_pohlig_hellman, &instance)
}

library_benchmark_group!(name = composite_order, benchmarks = [pohlig_hellman]);

main!(
    library_benchmark_groups = factorization,
    n_order_group,
    prime_order,
    composite_order
);
