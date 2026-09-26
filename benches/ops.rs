//! Instruction-count benchmarks (Valgrind/Callgrind) of the building blocks of `discrete_log`.
//!
//! Every input is deterministic (fixed seeds), and the setup work (generating instances,
//! computing the order given to the algorithms) is not measured. Results are checked: a wrong
//! answer makes the benchmark fail instead of reporting a fake speedup.
//!
//! Run with `cargo bench --features bench --bench ops` (requires Valgrind and `gungraun-runner`).

mod common;

use common::{DIGITS_108, Instance, SEED, check_factors, int, prime_bits};
use discrete_logarithm::{
    Error,
    bench::{element_order_with_factors, fast_factor, is_smooth},
    discrete_log_index_calculus_with_seed, discrete_log_pohlig_hellman,
    discrete_log_pollard_rho_with_seed, discrete_log_shanks_steps, n_order,
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

/// Small factor times a prime larger than every trial divisor: trial division stops at the
/// square root, the cofactor is recognized as prime.
fn large_prime_cofactor() -> Integer {
    prime_bits(40, SEED) * 587u32
}

/// Two primes above the trial division bound: only Pollard's rho splits them.
fn two_large_primes() -> Integer {
    Integer::from(33554467u32) * 33554473u32
}

/// A large prime squared: Pollard's rho never splits it, the square root is factored instead.
fn large_prime_squared() -> Integer {
    Integer::from(33554467u32).square()
}

/// `p - 1` of the 108-digit instance: prime powers just above the trial division bound, the
/// numbers the order of a base modulo a large prime is made of.
fn medium_prime_powers() -> Integer {
    int(DIGITS_108.0) - 1u32
}

#[library_benchmark]
#[bench::small(int("587"))]
#[bench::smooth(smooth())]
#[bench::large_prime_cofactor(large_prime_cofactor())]
#[bench::two_large_primes(two_large_primes())]
#[bench::large_prime_squared(large_prime_squared())]
#[bench::medium_prime_powers(medium_prime_powers())]
fn factor(n: Integer) -> HashMap<Integer, usize> {
    let factors = black_box(fast_factor(black_box(&n)));
    check_factors(&n, &factors);
    factors
}

/// Factor base of the index calculus for `n = 999231337607`: the primes below 508.
fn factorbase() -> Vec<u32> {
    Primes::all()
        .take_while(|&p| p < 508)
        .map(|p| p as u32)
        .collect()
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
fn index_calculus_smoothness(input: (Integer, Vec<u32>)) -> Option<Vec<u32>> {
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

/// `(n, b, factorization of n)`: what the order of `b` is computed from.
fn order_input(n: &str, b: &str) -> (Integer, Integer, HashMap<Integer, usize>) {
    let n = int(n);
    let n_factors = fast_factor(&n);
    (n, int(b), n_factors)
}

// The order of the base and its factorization, computed in one pass: what `discrete_log` spends
// most of its time on before it even chooses an algorithm.
#[library_benchmark]
#[bench::composite(order_input("32942478", "11"))]
// A modulus divisible by 8, where the group of units is not cyclic: the order divides the exact
// Carmichael lambda `2**38` of `2**40`, not `phi = 2**39`, and 3 is a unit modulo it.
#[bench::power_of_two(order_input("1099511627776", "3"))]
#[bench::digits_108(order_input(DIGITS_108.0, DIGITS_108.2))]
fn element_order(
    input: (Integer, Integer, HashMap<Integer, usize>),
) -> (Integer, HashMap<Integer, usize>) {
    let (n, b, n_factors) = black_box(&input);
    let (order, order_factors) = black_box(element_order_with_factors(n, b, n_factors));
    assert_eq!(
        b.clone().pow_mod(&order, n).unwrap(),
        1,
        "wrong order {order} of {b} modulo {n}"
    );
    check_factors(&order, &order_factors);
    (order, order_factors)
}

library_benchmark_group!(name = n_order_group, benchmarks = [order, element_order]);

type Algorithm = fn(&Integer, &Integer, &Integer, Option<&Integer>) -> Result<Integer, Error>;

/// Pollard's rho, seeded: the walks are the same in every run.
fn pollard_rho(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    discrete_log_pollard_rho_with_seed(n, a, b, order, SEED)
}

/// Index calculus, seeded: the relations are looked for in the same order in every run.
fn index_calculus_seeded(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    discrete_log_index_calculus_with_seed(n, a, b, order, SEED)
}

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
// fastest, which is what the thresholds of the algorithm selection (`SHANKS_STEPS_ORDER` and
// `INDEX_CALCULUS_SLACK` of `lib.rs`) are set from. The order of a safe prime of `bits` bits has
// `bits - 1` of them, so a boundary on the order falls one bit above the case that brackets it.

// Baby-step giant-step stops at 40 bits: the order of a safe prime of 41 bits is above `MAX_ORDER`,
// the size its table stops fitting in memory at, and it refuses such an order outright.
#[library_benchmark]
#[bench::bits_28(Instance::safe_prime(28, SEED))]
#[bench::bits_32(Instance::safe_prime(32, SEED))]
#[bench::bits_34(Instance::safe_prime(34, SEED))]
#[bench::bits_36(Instance::safe_prime(36, SEED))]
#[bench::bits_40(Instance::safe_prime(40, SEED))]
fn prime_order_shanks_steps(instance: Instance) -> Vec<Integer> {
    solve(discrete_log_shanks_steps, &instance)
}

#[library_benchmark]
#[bench::bits_28(Instance::safe_prime(28, SEED))]
#[bench::bits_32(Instance::safe_prime(32, SEED))]
#[bench::bits_34(Instance::safe_prime(34, SEED))]
#[bench::bits_36(Instance::safe_prime(36, SEED))]
#[bench::bits_40(Instance::safe_prime(40, SEED))]
#[bench::bits_42(Instance::safe_prime(42, SEED))]
#[bench::bits_46(Instance::safe_prime(46, SEED))]
fn prime_order_pollard_rho(instance: Instance) -> Vec<Integer> {
    solve(pollard_rho, &instance)
}

// Index calculus is the algorithm chosen from a modulus of about 45 bits on
// (`INDEX_CALCULUS_SLACK`), which the 42-bit and 46-bit cases bracket.
#[library_benchmark]
#[bench::bits_28(Instance::safe_prime(28, SEED))]
#[bench::bits_32(Instance::safe_prime(32, SEED))]
#[bench::bits_34(Instance::safe_prime(34, SEED))]
#[bench::bits_36(Instance::safe_prime(36, SEED))]
#[bench::bits_40(Instance::safe_prime(40, SEED))]
#[bench::bits_42(Instance::safe_prime(42, SEED))]
#[bench::bits_46(Instance::safe_prime(46, SEED))]
fn prime_order_index_calculus(instance: Instance) -> Vec<Integer> {
    solve(index_calculus_seeded, &instance)
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
    solve(index_calculus_seeded, &instance)
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
// The order of 2 modulo `1009**3` is `2**3 * 3**2 * 7 * 1009**2`: the two digits of `1009**2` are
// solved over one shared table of baby steps, where the digits of `2**3` and `3**2` fall below the
// exhaustive search threshold and are solved one by one.
#[bench::prime_power_order(Instance::new("1027243729", "2", &["177424694"]))]
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
