use std::collections::HashMap;

use rug::{Integer, ops::Pow};

use crate::{
    Algorithm, Error, algorithm_for, check_factorization, check_modulus,
    crt::crt,
    discrete_log_with_prime_order,
    factor::fast_factor,
    modular::{BigRing, ModRing, WordRing},
    n_order,
    shanks_steps::BabyStepTable,
};

/// Pohlig-Hellman algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// In order to compute the discrete logarithm, the algorithm takes advantage of the factorization of the group order. It is more efficient when the group order factors into many small primes.
///
/// A modulus that is not positive is refused with [`Error::InvalidModulus`], and an `order` that is
/// not positive with [`Error::InvalidOrder`]. Modulo 1 every residue is 0, so the logarithm is 0.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation. It
/// must be the order of `b` modulo `n` for the result to be right at all: a proper multiple of it
/// can make this return `Err(Error::LogDoesNotExist)` although a logarithm exists (`n = 7`,
/// `b = 2`, `a = 4`, `order = 9`), because the residues the prime powers of that multiple give have
/// no common solution. Sympy behaves the same way.
pub fn discrete_log_pohlig_hellman(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::new());
    }
    let order = match order {
        Some(order) => order.clone(),
        None => n_order(b, n)?,
    };
    // `fast_factor` of a number below 2 is the empty map, which is the factorization of 1 and of
    // nothing else: an order that is not positive has none, and is not an order of anything.
    if order < 1 {
        return Err(Error::InvalidOrder);
    }
    // The factorization is this crate's own: nothing to check about it.
    solve_with_factors(n, a, b, &order, &fast_factor(&order))
}

/// Pohlig-Hellman algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// Same as [`discrete_log_pohlig_hellman`] with the prime factorization of the order known: it
/// is the factorization the algorithm is built on, so nothing is left to compute.
///
/// `order_factors` must be the prime factorization of `order`, its primes mapped to their
/// exponents. Anything else is refused with [`Error::InvalidOrder`], the prime powers being
/// multiplied back together first: every projection, every digit and the final recombination is
/// built out of those prime powers, and a map belonging to another number describes no group.
pub fn discrete_log_pohlig_hellman_with_factors(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    order_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::new());
    }
    check_factorization(order, order_factors)?;
    solve_with_factors(n, a, b, order, order_factors)
}

/// [`discrete_log_pohlig_hellman_with_factors`] with `order_factors` already known to be the prime
/// factorization of `order`, which the callers inside the crate compute themselves.
pub(crate) fn solve_with_factors(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    order_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    debug_assert!(
        check_factorization(order, order_factors).is_ok(),
        "not the factorization of the order"
    );
    check_modulus(n)?;
    // Modulo 1 every residue is 0, so every exponent is a logarithm and 0 is the smallest one.
    if *n == 1 {
        return Ok(Integer::new());
    }

    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);

    // The word-size ring whenever the modulus fits in one, GMP otherwise.
    match WordRing::from_modulus(n) {
        Some(ring) => solve_in_ring(n, &ring, &a, &b, order, order_factors),
        None => {
            let ring = BigRing::new(n).expect("a positive modulus");
            solve_in_ring(n, &ring, &a, &b, order, order_factors)
        }
    }
}

/// [`solve`] in `ring`, `a` and `b` being already reduced modulo the modulus of it.
fn solve_in_ring<R: ModRing>(
    n: &Integer,
    ring: &R,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    order_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    let a = ring.from_integer(a);
    let b = ring.from_integer(b);

    // One residue per prime power of the order, from a discrete logarithm modulo each prime.
    let mut residues = Vec::with_capacity(order_factors.len());
    let mut modulli = Vec::with_capacity(order_factors.len());

    // `pi` is at least 2 and `ri` at least 1, checked by the caller: every prime power below is a
    // real one, and the indices into `powers` are all inside it.
    for (pi, ri) in order_factors {
        // `pi**0, ..., pi**(ri - 1)`: the weight of every digit, and the exponents the projections
        // below are taken to, computed once instead of once per digit.
        let powers: Vec<Integer> = (0..*ri as u32).map(|j| pi.clone().pow(j)).collect();
        let pe = Integer::from(&powers[*ri - 1] * pi);

        // The textbook projection into the subgroup of order `pi**ri`: from here on every exponent
        // is the size of that prime power instead of the size of the whole order.
        let cofactor = Integer::from(order / &pe);
        let a_p = ring.pow(&a, &cofactor);
        let b_p = ring.pow(&b, &cofactor);
        // `b_p**(pi**(ri - 1))`, of order `pi`: the base of every digit, and the same for all of
        // them where sympy rebuilds it from `b` for each one.
        let base = if *ri == 1 {
            b_p.clone()
        } else {
            ring.pow(&b_p, &powers[*ri - 1])
        };

        residues.push(digits(n, ring, pi, &powers, &a_p, &b_p, &base)?);
        modulli.push(pe);
    }

    crt(&residues, &modulli).ok_or(Error::LogDoesNotExist)
}

/// The residue of the logarithm modulo `pi**ri`, one digit of base `pi` at a time.
///
/// `a_p` and `b_p` are `a` and `b` projected into the subgroup of order `pi**ri`, `base` is
/// `b_p**(pi**(ri - 1))` and `powers` holds `pi**0, ..., pi**(ri - 1)`.
fn digits<R: ModRing>(
    n: &Integer,
    ring: &R,
    pi: &Integer,
    powers: &[Integer],
    a_p: &R::Elem,
    b_p: &R::Elem,
    base: &R::Elem,
) -> Result<Integer, Error> {
    // One table of baby steps for every digit, where sympy rebuilds a whole one per digit.
    let table = match shared_baby_steps(n, pi, powers.len()) {
        Some(order) => {
            Some(BabyStepTable::new(ring, base, order).ok_or(Error::NotRelativelyPrime)?)
        }
        None => None,
    };

    let mut residue = Integer::new();
    for (j, weight) in powers.iter().enumerate() {
        // `(a_p * b_p**-residue)**(pi**(ri - 1 - j))`, whose logarithm in base `base` is the digit
        // `j` of the residue: the digits already found are divided out first.
        let mut target = a_p.clone();
        if residue != 0 {
            let power = ring.pow(b_p, &residue);
            let inverse = ring.invert(&power).ok_or(Error::NotRelativelyPrime)?;
            ring.mul_assign(&mut target, &inverse);
        }
        // The last digit is already in the subgroup of order `pi`: its exponent is `pi**0`.
        if j + 1 < powers.len() {
            ring.pow_assign(&mut target, &powers[powers.len() - 1 - j]);
        }

        let digit = match &table {
            Some(table) => table
                .log(&target)
                .map(Integer::from)
                .ok_or(Error::LogDoesNotExist)?,
            // `pi` is prime: the sub-problem never comes back to Pohlig-Hellman.
            None => discrete_log_with_prime_order(
                n,
                &ring.to_integer(&target),
                &ring.to_integer(base),
                pi,
            )?,
        };
        residue += digit * weight;
    }

    Ok(residue)
}

/// The order of the sub-problems of `pi**ri` when their digits are worth a shared table of baby
/// steps, `None` when they are not.
///
/// The table stands in for the baby-step giant-step [`discrete_log_with_prime_order`] would run
/// for each digit, and for nothing else: a single digit shares nothing, and a prime the dispatcher
/// would solve any other way is left to it. [`algorithm_for`] is asked which that is, so that the
/// algorithm selection cannot drift apart from the one the dispatcher makes.
fn shared_baby_steps(n: &Integer, pi: &Integer, ri: usize) -> Option<u64> {
    if ri < 2 {
        return None;
    }
    // Every order baby-step giant-step is chosen for is below `MAX_ORDER`, so it fits in a `u64`.
    let order = pi.to_u64()?;
    match algorithm_for(n, pi) {
        Algorithm::ShanksSteps => Some(order),
        Algorithm::TrialMul | Algorithm::IndexCalculus | Algorithm::PollardRho => None,
    }
}

#[cfg(test)]
mod tests {
    use rug::ops::Pow;

    use super::*;

    #[test]
    fn pohlig_hellman() {
        assert_eq!(
            discrete_log_pohlig_hellman(
                &98376431.into(),
                &(Integer::from(11).pow(9)),
                &11.into(),
                None
            )
            .unwrap(),
            9
        );
        assert_eq!(
            discrete_log_pohlig_hellman(
                &78723213.into(),
                &(Integer::from(11).pow(31)),
                &11.into(),
                None
            )
            .unwrap(),
            31
        );
        assert_eq!(
            discrete_log_pohlig_hellman(
                &32942478.into(),
                &(Integer::from(11).pow(98)),
                &11.into(),
                None
            )
            .unwrap(),
            98
        );
        assert_eq!(
            discrete_log_pohlig_hellman(
                &14789363.into(),
                &(Integer::from(11).pow(444)),
                &11.into(),
                None
            )
            .unwrap(),
            444
        );
    }

    #[test]
    fn known_order_factors() {
        // The factorization of the order is the only thing the algorithm needs.
        let n = Integer::from(32942478);
        let order = Integer::from(2745206);
        let order_factors = fast_factor(&order);
        assert_eq!(
            discrete_log_pohlig_hellman_with_factors(
                &n,
                &(Integer::from(11).pow(98)),
                &11.into(),
                &order,
                &order_factors
            )
            .unwrap(),
            98
        );
        assert_eq!(
            discrete_log_pohlig_hellman(&n, &(Integer::from(11).pow(98)), &11.into(), Some(&order))
                .unwrap(),
            98
        );
    }

    #[test]
    fn repeated_prime_of_the_order() {
        // The order of 2 modulo `1009**3` is `2**3 * 3**2 * 7 * 1009**2`: the two digits of
        // `1009**2` share one table of baby steps (1009 is above the exhaustive search
        // threshold), where the digits of `2**3` and `3**2` are still solved one by one.
        let n = Integer::from(1009).pow(3u32);
        let order = n_order(&2.into(), &n).unwrap();
        let order_factors = fast_factor(&order);
        assert_eq!(order_factors.get(&Integer::from(1009)), Some(&2));

        for x in [0u32, 1, 2, 1009, 1_018_081, 123_456_789] {
            let a = Integer::from(2).pow_mod(&Integer::from(x), &n).unwrap();
            assert_eq!(
                discrete_log_pohlig_hellman_with_factors(&n, &a, &2.into(), &order, &order_factors)
                    .unwrap(),
                x,
                "log of 2**{x} in base 2 modulo {n}"
            );
        }

        // 11 is not a power of 2 modulo `n`: the shared table must not answer a logarithm that
        // does not exist.
        assert_eq!(
            discrete_log_pohlig_hellman_with_factors(
                &n,
                &11.into(),
                &2.into(),
                &order,
                &order_factors
            ),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn base_that_is_not_invertible() {
        // The order of `2**38` is reached before the digits of it are exhausted, and dividing the
        // digits found so far out of the target needs an inverse of the base that does not exist.
        let n = Integer::from(2).pow(40u32) * 3u32;
        let order = Integer::from(2).pow(38u32) * 3u32;
        let order_factors = HashMap::from([(Integer::from(2), 38), (Integer::from(3), 1)]);
        assert_eq!(
            discrete_log_pohlig_hellman_with_factors(
                &n,
                &4.into(),
                &2.into(),
                &order,
                &order_factors
            ),
            Err(Error::NotRelativelyPrime)
        );
    }

    #[test]
    fn order_factors_validation() {
        // Every one of these used to panic: a division by zero, a negative exponent handed to GMP,
        // or a `debug_assert` of the ring arithmetic.
        let big = Integer::from(1) << 64u32;
        let big_even = (Integer::from(1) << 70u32) + 4u32;
        for (n, order, factors) in [
            // A prime of 0, which nothing can be divided by.
            (Integer::from(7), Integer::from(6), vec![(0, 1)]),
            // A prime of 1, whose powers are all 1.
            (Integer::from(7), Integer::from(6), vec![(1, 1)]),
            // Negative primes, one of them with a modulus GMP holds.
            (Integer::from(7), Integer::from(6), vec![(-3, 1)]),
            (big_even.clone(), Integer::from(6), vec![(-3, 1)]),
            (big_even, Integer::from(6), vec![(2, 1), (-3, 1)]),
            // Orders that are not positive.
            (big.clone(), Integer::from(-5), vec![(2, 1)]),
            (big.clone(), Integer::from(-2), vec![(2, 1)]),
            (big.clone(), Integer::from(-1), vec![(1, 1)]),
            (big.clone(), Integer::from(0), vec![(2, 1)]),
            // An exponent of 0, and one no prime power can have.
            (big.clone(), Integer::from(1), vec![(2, 0)]),
            (big.clone(), Integer::from(6), vec![(2, usize::MAX)]),
            // Factorizations of another number: the prime powers do not multiply back to the order.
            (big.clone(), Integer::from(6), vec![(2, 2)]),
            (big.clone(), Integer::from(6), vec![(3, 1)]),
            (big, Integer::from(6), vec![(2, 1), (3, 1), (5, 1)]),
        ] {
            let order_factors: HashMap<Integer, usize> = factors
                .iter()
                .map(|&(p, e)| (Integer::from(p), e))
                .collect();
            assert_eq!(
                discrete_log_pohlig_hellman_with_factors(
                    &n,
                    &4.into(),
                    &2.into(),
                    &order,
                    &order_factors
                ),
                Err(Error::InvalidOrder),
                "n = {n}, order = {order}, factors = {factors:?}"
            );
        }

        // And the factorization that does belong to the order still solves.
        assert_eq!(
            discrete_log_pohlig_hellman_with_factors(
                &7.into(),
                &4.into(),
                &2.into(),
                &3.into(),
                &HashMap::from([(Integer::from(3), 1)])
            )
            .unwrap(),
            2
        );

        // An order that is not positive is refused through the other entry point too, which factors
        // it itself: there is no factorization of it to give.
        for order in [-5, -1, 0] {
            assert_eq!(
                discrete_log_pohlig_hellman(&7.into(), &4.into(), &2.into(), Some(&order.into())),
                Err(Error::InvalidOrder),
                "order {order}"
            );
        }
    }

    #[test]
    fn order_that_is_a_multiple_of_the_real_one() {
        // The order of 2 modulo 7 is 3, and 9 is not a multiple of it: the residues of the prime
        // powers of 9 have no common solution, so a logarithm that exists is not found. Sympy
        // answers the same way, and this is what the documented precondition on `order` is about.
        assert_eq!(
            discrete_log_pohlig_hellman(&7.into(), &4.into(), &2.into(), Some(&9.into())),
            Err(Error::LogDoesNotExist)
        );
        // A real multiple of it still gives a solution, and not always the smallest one.
        let x = discrete_log_pohlig_hellman(&7.into(), &4.into(), &2.into(), Some(&6.into()))
            .unwrap()
            .to_u32()
            .unwrap();
        assert_eq!(
            Integer::from(2)
                .pow_mod(&Integer::from(x), &7.into())
                .unwrap(),
            4
        );
    }

    #[test]
    fn modulus_validation() {
        // Modulo 1 every number is 0.
        assert_eq!(
            discrete_log_pohlig_hellman(&1.into(), &0.into(), &0.into(), None).unwrap(),
            0
        );
        for n in [0, -1, -4] {
            assert_eq!(
                discrete_log_pohlig_hellman(&n.into(), &1.into(), &3.into(), Some(&10.into())),
                Err(Error::InvalidModulus)
            );
        }
    }
}
