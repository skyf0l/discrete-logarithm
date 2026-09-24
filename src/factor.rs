use std::{collections::HashMap, sync::OnceLock};

use primal::Primes;
use rug::{integer::IsPrime, ops::Pow, Assign, Integer};

/// Trial division is done by every prime below this bound.
const TRIAL_DIVISION_BOUND: u32 = 1 << 16;

/// Number of Miller-Rabin rounds used to tell primes from composites.
const PRIMALITY_REPS: u32 = 30;

/// Primes below [`TRIAL_DIVISION_BOUND`], sieved once for the whole process.
fn trial_primes() -> &'static [u32] {
    static PRIMES: OnceLock<Vec<u32>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        Primes::all()
            .take_while(|&p| p < TRIAL_DIVISION_BOUND as usize)
            .map(|p| p as u32)
            .collect()
    })
}

/// Returns the prime factorization of `n`, as a map of primes to their exponents.
///
/// `n` is trial divided by the primes below 65536, then the remaining cofactor is split with
/// Brent's variant of Pollard's rho. Numbers with two large prime factors (both above ~10^12)
/// are out of reach: factoring them is as hard as the general factoring problem.
///
/// Numbers smaller than 2 have no prime factors: the map is empty.
pub fn fast_factor(n: &Integer) -> HashMap<Integer, usize> {
    let mut factors = HashMap::new();
    if *n <= 1 {
        return factors;
    }

    let mut n = n.clone();
    for &p in trial_primes() {
        // Every prime below `p` has been divided out: a cofactor below `p**2` is prime.
        if n < u64::from(p) * u64::from(p) {
            break;
        }
        if n.is_divisible_u(p) {
            let mut exponent = 0;
            while n.is_divisible_u(p) {
                n.div_exact_u_mut(p);
                exponent += 1;
            }
            factors.insert(Integer::from(p), exponent);
        }
    }

    factor_cofactor(n, &mut factors);
    factors
}

/// Adds the prime factorization of `n` to `factors`, `n` having no prime factor below
/// [`TRIAL_DIVISION_BOUND`].
fn factor_cofactor(n: Integer, factors: &mut HashMap<Integer, usize>) {
    if n == 1 {
        return;
    }
    if n.is_probably_prime(PRIMALITY_REPS) != IsPrime::No {
        *factors.entry(n).or_insert(0) += 1;
        return;
    }

    // Pollard's rho is slow on `p**k`, and the roots of a power are much easier to factor.
    if let Some((root, exponent)) = perfect_power(&n) {
        let mut root_factors = HashMap::new();
        factor_cofactor(root, &mut root_factors);
        for (p, e) in root_factors {
            *factors.entry(p).or_insert(0) += e * exponent;
        }
        return;
    }

    let divisor = pollard_rho(&n);
    let cofactor = Integer::from(&n / &divisor);
    factor_cofactor(divisor, factors);
    factor_cofactor(cofactor, factors);
}

/// Returns `(root, exponent)` with `root ** exponent == n` and `exponent` prime, if `n` is a
/// perfect power.
fn perfect_power(n: &Integer) -> Option<(Integer, usize)> {
    if !n.is_perfect_power() {
        return None;
    }
    // A perfect power is a perfect `p`-th power for at least one prime `p`.
    for &p in trial_primes() {
        if u64::from(p) > u64::from(n.significant_bits()) {
            break;
        }
        let root = Integer::from(n.root_ref(p));
        if root.clone().pow(p) == *n {
            return Some((root, p as usize));
        }
    }
    None
}

/// Returns a non trivial divisor of the composite `n`, with Brent's variant of Pollard's rho.
///
/// `n` must have at least two distinct prime factors, else the cycles never split it.
fn pollard_rho(n: &Integer) -> Integer {
    // The polynomials `x**2 + c` are tried in order, so the same number always splits the same
    // way: no randomness, reproducible running times.
    let mut c = 1;
    loop {
        if let Some(divisor) = pollard_rho_cycle(n, c) {
            return divisor;
        }
        c += 1;
    }
}

/// One attempt at splitting `n` with the cycles of `x -> x**2 + c`, `None` on failure.
fn pollard_rho_cycle(n: &Integer, c: u32) -> Option<Integer> {
    // Products of that many differences share a single gcd: gcd is much slower than a
    // multiplication modulo `n`.
    const BATCH: u64 = 128;

    let mut y = Integer::from(2);
    let mut x = Integer::new();
    let mut last = Integer::new();
    let mut diff = Integer::new();
    let mut product = Integer::from(1);
    let mut divisor = Integer::from(1);
    // Length of the cycle looked for, doubled at each round.
    let mut range = 1u64;

    while divisor == 1 {
        x.assign(&y);
        for _ in 0..range {
            rho_step(&mut y, c, n);
        }
        let mut done = 0;
        while done < range && divisor == 1 {
            last.assign(&y);
            for _ in 0..BATCH.min(range - done) {
                rho_step(&mut y, c, n);
                diff.assign(&x - &y);
                diff.abs_mut();
                product *= &diff;
                product %= n;
            }
            divisor.assign(&product);
            divisor.gcd_mut(n);
            done += BATCH;
        }
        range *= 2;
    }

    if divisor == *n {
        // The batch multiplied the factor by its cofactor: walk the last batch step by step.
        loop {
            rho_step(&mut last, c, n);
            diff.assign(&x - &last);
            diff.abs_mut();
            divisor.assign(&diff);
            divisor.gcd_mut(n);
            if divisor != 1 {
                break;
            }
        }
    }

    if divisor == *n {
        None
    } else {
        Some(divisor)
    }
}

/// `y = y**2 + c (mod n)`.
fn rho_step(y: &mut Integer, c: u32, n: &Integer) {
    y.square_mut();
    *y += c;
    *y %= n;
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    /// Factors of `n`, sorted by prime.
    fn factor(n: &str) -> Vec<(Integer, usize)> {
        let mut factors = fast_factor(&Integer::from_str(n).unwrap())
            .into_iter()
            .collect::<Vec<_>>();
        factors.sort();
        factors
    }

    fn check(n: &str, expected: &[(&str, usize)]) {
        let expected = expected
            .iter()
            .map(|(p, e)| (Integer::from_str(p).unwrap(), *e))
            .collect::<Vec<_>>();
        assert_eq!(factor(n), expected, "wrong factorization of {n}");
    }

    #[test]
    fn factorization() {
        // Numbers without prime factors.
        check("0", &[]);
        check("1", &[]);
        for n in ["-1", "-2", "-100"] {
            check(n, &[]);
        }

        // Trial division.
        check("2", &[("2", 1)]);
        check("4", &[("2", 2)]);
        check("587", &[("587", 1)]);
        check("1000000", &[("2", 6), ("5", 6)]);
        check("32942478", &[("2", 1), ("3", 1), ("5490413", 1)]);

        // Cofactor larger than every trial divisor.
        check("65537", &[("65537", 1)]);
        check("4295491591", &[("65537", 1), ("65543", 1)]);
    }

    #[test]
    fn factorization_large_factors() {
        // Both factors above the trial division bound: Pollard's rho.
        check("1125902456980891", &[("33554467", 1), ("33554473", 1)]);
        // A large prime squared: Pollard's rho never splits it, the root is factored instead.
        check("1125902255654089", &[("33554467", 2)]);
        // A composite squared.
        check(
            "1267656342635607108898739153881",
            &[("33554467", 2), ("33554473", 2)],
        );
        // Three large prime factors.
        check(
            "37779097370472678002173",
            &[("33554467", 1), ("33554473", 1), ("33554503", 1)],
        );
    }

    #[test]
    fn factorization_products_are_exact() {
        for n in [
            "982451653",
            "99999999999999999999999999",
            "170141183460469231731687303715884105727",
            "618970019642690137449562111",
        ] {
            let n = Integer::from_str(n).unwrap();
            let factors = fast_factor(&n);
            let product = factors.iter().fold(Integer::from(1), |acc, (p, e)| {
                acc * p.clone().pow(*e as u32)
            });
            assert_eq!(product, n, "wrong factorization of {n}: {factors:?}");
            for p in factors.keys() {
                assert_ne!(p.is_probably_prime(PRIMALITY_REPS), IsPrime::No);
            }
        }
    }
}
