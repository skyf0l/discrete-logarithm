use std::{collections::HashMap, iter, sync::OnceLock};

use rug::{Assign, Integer, integer::IsPrime, ops::Pow};

use crate::modular::{ModRing, WordRing};

/// Trial division is done by every prime below this bound.
const TRIAL_DIVISION_BOUND: u32 = 1 << 16;

/// Number of Miller-Rabin rounds used to tell primes from composites, here and in the algorithm
/// selection of [`crate::discrete_log`].
///
/// GMP runs a Baillie-PSW test first, which no composite is known to pass, then this many rounds
/// with random bases, each of which a composite passes with probability at most 1/4: 25 is the
/// usual choice, and the few modular powerings the last rounds add cost nothing next to a
/// factorization or a discrete logarithm.
pub(crate) const PRIMALITY_REPS: u32 = 30;

/// The primes below 256.
///
/// Trial division stops as soon as the cofactor is below the square of the next divisor, so these
/// are the only ones a number below `256**2` can need: the rest of the table is sieved only for
/// the larger ones.
const SMALL_PRIMES: [u32; 54] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
];

/// Number of primes between 256 and [`TRIAL_DIVISION_BOUND`]: `pi(65536) = 6542`, less the 54 of
/// [`SMALL_PRIMES`].
const SIEVED_PRIMES: usize = 6488;

/// The primes from 256 to [`TRIAL_DIVISION_BOUND`], sieved once for the whole process.
fn sieved_primes() -> &'static [u32] {
    static PRIMES: OnceLock<Vec<u32>> = OnceLock::new();
    PRIMES.get_or_init(|| {
        // Sieve of Eratosthenes over the odd numbers alone: `index` stands for `2 * index + 1`.
        let odds = TRIAL_DIVISION_BOUND as usize / 2;
        let mut composite = vec![false; odds];
        // Sized once for all the primes it will hold, the count of them below the bound being known:
        // growing it instead costs a capacity check per push, which shows up as noise on every
        // benchmark that factors anything. The assertion below is what keeps the two in step.
        let mut primes = Vec::with_capacity(SIEVED_PRIMES);
        // Index 0 stands for 1, which is neither prime nor a useful divisor.
        for index in 1..odds {
            if composite[index] {
                continue;
            }
            let prime = 2 * index + 1;
            if prime > SMALL_PRIMES[SMALL_PRIMES.len() - 1] as usize {
                primes.push(prime as u32);
            }
            // Every odd multiple of `prime` below its square has a smaller prime factor.
            let mut multiple = prime * prime / 2;
            while multiple < odds {
                composite[multiple] = true;
                multiple += prime;
            }
        }
        debug_assert_eq!(
            primes.len(),
            SIEVED_PRIMES,
            "the sieve holds another count of primes"
        );
        primes
    })
}

/// Every trial divisor, in increasing order: the primes below [`TRIAL_DIVISION_BOUND`].
///
/// The sieve of the primes above 256 is built on first use, so the numbers split by
/// [`SMALL_PRIMES`] alone never pay for it.
fn trial_primes() -> impl Iterator<Item = u32> {
    SMALL_PRIMES
        .into_iter()
        .chain(iter::once_with(sieved_primes).flatten().copied())
}

/// Returns the prime factorization of `n`, as a map of primes to their exponents.
///
/// `n` is trial divided by the primes below 65536, then the remaining cofactor is split with
/// Brent's variant of Pollard's rho. Two prime factors of about 2^50 each take a couple of
/// seconds, 2^56 each about half a minute, and beyond that factoring `n` is the general factoring
/// problem and out of reach.
///
/// Numbers smaller than 2 have no prime factors: the map is empty.
pub fn fast_factor(n: &Integer) -> HashMap<Integer, usize> {
    let mut factors = HashMap::new();
    if *n <= 1 {
        return factors;
    }

    let mut n = n.clone();
    match n.to_u64() {
        // A remainder of a `u64` is a single instruction, where GMP pays a call and a size check
        // for each of the thousands of trial divisors.
        Some(word) => n.assign(trial_divide_word(word, &mut factors)),
        None => trial_divide(&mut n, &mut factors),
    }

    factor_cofactor(n, &mut factors);
    factors
}

/// Divides the trial divisors out of `n`, adding those that divide it to `factors`.
fn trial_divide(n: &mut Integer, factors: &mut HashMap<Integer, usize>) {
    for p in trial_primes() {
        // Every prime below `p` has been divided out: a cofactor below `p**2` is prime.
        if *n < u64::from(p) * u64::from(p) {
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
}

/// Divides the trial divisors out of the `u64` `n`, adding those that divide it to `factors`, and
/// returns the cofactor.
fn trial_divide_word(mut n: u64, factors: &mut HashMap<Integer, usize>) -> u64 {
    for p in trial_primes() {
        let p = u64::from(p);
        // Every prime below `p` has been divided out: a cofactor below `p**2` is prime.
        if n < p * p {
            break;
        }
        if n % p == 0 {
            let mut exponent = 0;
            while n % p == 0 {
                n /= p;
                exponent += 1;
            }
            factors.insert(Integer::from(p), exponent);
        }
    }
    n
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
    for p in trial_primes() {
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

/// Products of that many differences share a single gcd: a gcd is much slower than a
/// multiplication modulo `n`.
const RHO_BATCH: u64 = 128;

/// One attempt at splitting `n` with the cycles of `x -> x**2 + c`, `None` on failure.
fn pollard_rho_cycle(n: &Integer, c: u32) -> Option<Integer> {
    // The cofactors left by trial division nearly always fit in a `u64`, where a whole walk runs
    // without a single allocation.
    match n.to_u64() {
        Some(word) => pollard_rho_cycle_word(word, c).map(Integer::from),
        None => pollard_rho_cycle_big(n, c),
    }
}

/// One attempt at splitting the `u64` `n` with the cycles of `x -> x**2 + c`, `None` on failure.
fn pollard_rho_cycle_word(n: u64, c: u32) -> Option<u64> {
    let ring = WordRing::new(n).expect("a composite is not zero");
    let c = ring.from_u64(u64::from(c));

    let mut y = ring.from_u64(2);
    // The point the walk is compared to, and the point the current batch started at.
    let mut x = y;
    let mut last = y;
    let mut product = ring.one();
    let mut divisor = 1;
    // Length of the cycle looked for, doubled at each round.
    let mut range = 1u64;

    while divisor == 1 {
        x = y;
        for _ in 0..range {
            y = rho_step_word(&ring, y, c);
        }
        let mut done = 0;
        while done < range && divisor == 1 {
            last = y;
            for _ in 0..RHO_BATCH.min(range - done) {
                y = rho_step_word(&ring, y, c);
                // A difference and its opposite have the same gcd with `n`: no need to take the
                // absolute value of `x - y`.
                product = ring.mul(&product, &ring.sub(&x, &y));
            }
            divisor = gcd(ring.to_u64(product), n);
            done += RHO_BATCH;
        }
        // Saturating: 63 doublings wrap a `u64` to zero, which turns the loop into a walk of no
        // steps that never ends. No factorization gets that far, and a walk of `u64::MAX` steps is
        // the honest meaning of a range that cannot grow any more.
        range = range.saturating_mul(2);
    }

    if divisor == n {
        // The batch multiplied the factor by its cofactor: walk the last batch step by step.
        loop {
            last = rho_step_word(&ring, last, c);
            divisor = gcd(ring.to_u64(ring.sub(&x, &last)), n);
            if divisor != 1 {
                break;
            }
        }
    }

    (divisor != n).then_some(divisor)
}

/// One attempt at splitting `n` with the cycles of `x -> x**2 + c`, `None` on failure.
fn pollard_rho_cycle_big(n: &Integer, c: u32) -> Option<Integer> {
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
            for _ in 0..RHO_BATCH.min(range - done) {
                rho_step(&mut y, c, n);
                diff.assign(&x - &y);
                diff.abs_mut();
                product *= &diff;
                product %= n;
            }
            divisor.assign(&product);
            divisor.gcd_mut(n);
            done += RHO_BATCH;
        }
        // Saturating: 63 doublings wrap a `u64` to zero, which turns the loop into a walk of no
        // steps that never ends. No factorization gets that far, and a walk of `u64::MAX` steps is
        // the honest meaning of a range that cannot grow any more.
        range = range.saturating_mul(2);
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

    if divisor == *n { None } else { Some(divisor) }
}

/// `y = y**2 + c (mod n)`.
fn rho_step(y: &mut Integer, c: u32, n: &Integer) {
    y.square_mut();
    *y += c;
    *y %= n;
}

/// `y**2 + c` in `ring`, `c` being already an element of it.
#[inline]
fn rho_step_word(ring: &WordRing, y: u64, c: u64) -> u64 {
    ring.add(&ring.square(&y), &c)
}

/// The greatest common divisor of `a` and `b`.
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a
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

    #[test]
    fn prime_table() {
        // The table is exactly the primes below the bound, in increasing order.
        let expected: Vec<u32> = (2..TRIAL_DIVISION_BOUND)
            .filter(|&n| (2..n).take_while(|d| d * d <= n).all(|d| n % d != 0))
            .collect();
        assert_eq!(trial_primes().collect::<Vec<_>>(), expected);
        // Splitting the table in two must not lose or repeat a prime.
        assert_eq!(SMALL_PRIMES.last(), Some(&251));
        assert_eq!(sieved_primes().first(), Some(&257));
    }

    #[test]
    fn word_and_gmp_rho_agree() {
        // The two implementations of a cycle walk the same points and split a number the same
        // way, so the factorizations do not depend on the size of the cofactor.
        for n in [
            "4295491591",
            "1125902456980891",
            // Many small factors, and an even number: the ring has no Montgomery form for it.
            "18446744073709551615",
            "8590983182",
            // Two primes just below `2**32`, the largest product a word holds.
            "18446743979220271189",
        ] {
            let n = Integer::from_str(n).unwrap();
            for c in 1..4 {
                assert_eq!(
                    pollard_rho_cycle_word(n.to_u64().unwrap(), c).map(Integer::from),
                    pollard_rho_cycle_big(&n, c),
                    "n = {n}, c = {c}"
                );
            }
        }
    }

    #[test]
    fn word_gcd() {
        assert_eq!(gcd(0, 0), 0);
        assert_eq!(gcd(0, 7), 7);
        assert_eq!(gcd(7, 0), 7);
        assert_eq!(gcd(12, 18), 6);
        assert_eq!(gcd(u64::MAX, u64::MAX - 1), 1);
    }
}
