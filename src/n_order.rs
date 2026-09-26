use std::collections::HashMap;

use rug::{Integer, ops::Pow};

use crate::{Error, check_factorization, check_modulus, factor::fast_factor};

/// Returns the order of `a` modulo `n`.
///
/// The order of `a` modulo `n` is the smallest integer `k` such that `a**k` leaves a remainder of 1 with `n`.
///
/// Modulo 1 every number is 1, so the order is 1.
pub fn n_order(a: &Integer, n: &Integer) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::from(1));
    }
    // The factorization is this crate's own: nothing to check about it.
    order_from_factors(a, n, &fast_factor(n))
}

/// Returns the order of `a` modulo `n`.
///
/// The order of `a` modulo `n` is the smallest integer `k` such that `a**k` leaves a remainder of 1 with `n`.
///
/// If the prime factorization of `n` is known, it can be passed as `n_factors` to speed up the computation.
///
/// `n_factors` must be the prime factorization of `n`, its primes mapped to their exponents.
/// Anything else is refused with [`Error::InvalidOrder`]: the order is built out of that map one
/// prime power at a time, and a map belonging to another number describes no group this `a` lives
/// in. The primes themselves are taken on trust, the product check being all that is paid for; a
/// composite one that it lets through is refused where it is met instead.
pub fn n_order_with_factors(
    a: &Integer,
    n: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::from(1));
    }
    check_factorization(n, n_factors)?;
    order_from_factors(a, n, n_factors)
}

/// [`n_order_with_factors`] with `n_factors` already known to be the factorization of `n`.
fn order_from_factors(
    a: &Integer,
    n: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    let a = a.clone().modulo(n);
    // Trivial case, and the only one where the factorization of `n` is not needed.
    if a == 1 {
        return Ok(Integer::from(1));
    }
    if a.clone().gcd(n) != 1 {
        return Err(Error::NotRelativelyPrime);
    }

    // The order modulo `n` is the lcm of the orders modulo each prime power of `n`: the
    // exponents stay as small as the prime power they belong to.
    let mut order = Integer::from(1);
    for (p, e) in n_factors {
        order.lcm_mut(&prime_power_order(&a, p, *e)?);
    }
    Ok(order)
}

/// Returns the order of `a` modulo `p**e`, `a` and `p` being relatively prime.
///
/// `Err(Error::InvalidOrder)` when they are not, which only a `p` that divides no part of the
/// modulus can bring: the caller checked `gcd(a, n) == 1`, and that says nothing about a `p` taken
/// from the factorization of another number. The powers of such an `a` never reach 1, and the loop
/// below would multiply the order by `px` for as long as memory allowed.
fn prime_power_order(a: &Integer, p: &Integer, e: usize) -> Result<Integer, Error> {
    let pe = p.clone().pow(e as u32);
    // The units modulo `p**e` form a group of order `(p - 1) * p**(e - 1)`, so the order of `a`
    // divides it.
    let group_order = Integer::from(p - 1u32) * p.clone().pow(e as u32 - 1);
    let mut group_order_factors = fast_factor(&Integer::from(p - 1u32));
    if e > 1 {
        group_order_factors.insert(p.clone(), e - 1);
    }

    let mut order = Integer::from(1);
    for (px, ex) in &group_order_factors {
        // Remove `px**ex` from the exponent, then put back the powers of `px` needed to reach 1.
        let exponent = &group_order / px.clone().pow(*ex as u32);
        let mut x = a.clone().pow_mod(&exponent, &pe).unwrap();
        // The order of `x` divides `px**ex` whenever `a` is a unit modulo `p**e`, so `ex` of these
        // powerings reach 1: one more says that it is not one, and nothing here is an order.
        let mut added = 0;
        while x != 1 {
            if added == *ex {
                return Err(Error::InvalidOrder);
            }
            x = x.pow_mod(px, &pe).unwrap();
            order *= px;
            added += 1;
        }
    }
    Ok(order)
}

/// Returns the order of `b` modulo `n` and the prime factorization of that order.
///
/// A multiple of the order of `b` and its factorization are computed in one pass, then the primes
/// the order of `b` does not need are divided out: the factorization of the order comes for free,
/// and the primality test of [`n_order`] is not needed.
///
/// `b` must be reduced modulo `n`. When `b` and `n` are not relatively prime, the powers of `b`
/// are not a subgroup of the units and nothing can be divided out: the returned order is the
/// order of the group, which the callers only use as a bound.
///
/// `n_factors` must be the prime factorization of `n`: the public entry points check it before they
/// get here.
pub fn element_order_with_factors(
    n: &Integer,
    b: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> (Integer, HashMap<Integer, usize>) {
    // The order of a unit divides Carmichael's lambda, the lcm of the exponents of the units modulo
    // each prime power of `n`, which is often much smaller than their product `totient(n)` (5490412
    // instead of 10980824 for `n = 32942478`): the loop below starts from a smaller multiple of the
    // order and strips it in fewer modular exponentiations.
    //
    // Lambda only bounds the order of a unit though. The powers of a base sharing a factor with
    // the modulus are not a subgroup of the units, and the callers use the returned value as a
    // bound of an exhaustive search: the logarithm of 4 in base 2 modulo 12 is 2, where
    // `lambda(12) = 2` would stop that search before reaching it. The totient is what is kept for
    // that case, as sympy does.
    //
    // It is not a bound on how many distinct powers such a base has: the powers of 0 modulo 2 are
    // 1 and 0, two of them for a totient of 1, and 534 modulo 586 has `totient(586) + 1` of them.
    // A logarithm that only the last of those powers gives is therefore missed here, exactly as it
    // is in sympy, and the differential test holds the two together.
    let unit = Integer::from(b.gcd_ref(n)) == 1;
    // A modulus with a single prime power has no lcm to take: `px` divides no factor of `px - 1`,
    // so each prime appears once and the lcm is the product.
    let lcm = unit && n_factors.len() > 1;

    let mut factors: HashMap<Integer, usize> = HashMap::new();
    for (px, kx) in n_factors {
        // `px` divides no `py - 1` of the same prime power, so this is its whole contribution to
        // the order or to the exponent of the units modulo `px**kx`.
        let exponent = if unit && *px == 2 {
            two_power_exponent(*kx)
        } else {
            kx - 1
        };
        if exponent > 0 {
            keep_exponent(&mut factors, px.clone(), exponent, lcm);
        }
        for (py, ky) in fast_factor(&Integer::from(px - 1u32)) {
            keep_exponent(&mut factors, py, ky, lcm);
        }
    }
    let mut order = factors.iter().fold(Integer::from(1), |order, (p, e)| {
        order * p.clone().pow(*e as u32)
    });

    // The order of `b` divides the order of the group.
    let mut order_factors = HashMap::new();
    for (p, e) in &factors {
        let mut removed = 0;
        for _ in 0..*e {
            let smaller = Integer::from(&order / p);
            if b.clone().pow_mod(&smaller, n).unwrap() != 1 {
                break;
            }
            order = smaller;
            removed += 1;
        }
        if removed < *e {
            order_factors.insert(p.clone(), e - removed);
        }
    }
    (order, order_factors)
}

/// The exponent of 2 in the exponent of the group of units modulo `2**e`.
///
/// That group is the one place where the exponent of the units is not the order of the group: it is
/// trivial modulo 2, of order 2 modulo 4, and `Z/2 x Z/2**(e - 2)` above, so its exponent is
/// `2**(e - 2)` from `e = 3` on where the order is `2**(e - 1)`. Taking the order instead starts
/// every modulus divisible by 8 from a bound twice too large, and costs one modular exponentiation
/// per factor of 2 to strip it back down.
fn two_power_exponent(e: usize) -> usize {
    debug_assert!(e > 0, "a prime power of exponent 0");
    match e {
        1 => 0,
        2 => 1,
        _ => e - 2,
    }
}

/// Merges `p**e` into `factors`, an `lcm` keeping the largest exponent where a product adds them.
fn keep_exponent(factors: &mut HashMap<Integer, usize>, p: Integer, e: usize, lcm: bool) {
    let exponent = factors.entry(p).or_insert(0);
    *exponent = if lcm {
        (*exponent).max(e)
    } else {
        *exponent + e
    };
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    fn int(s: &str) -> Integer {
        Integer::from_str(s).unwrap()
    }

    #[test]
    fn orders() {
        assert_eq!(n_order(&2.into(), &13.into()).unwrap(), 12);
        for (a, res) in (1..=6).zip([1, 3, 6, 3, 6, 2]) {
            assert_eq!(n_order(&a.into(), &7.into()).unwrap(), res);
        }
        assert_eq!(n_order(&5.into(), &17.into()).unwrap(), 16);
        assert_eq!(
            n_order(&17.into(), &11.into()),
            n_order(&6.into(), &11.into())
        );
        assert_eq!(n_order(&101.into(), &119.into()).unwrap(), 6);
        assert_eq!(n_order(&3.into(), &2456747.into()).unwrap(), 1228373);
        assert_eq!(n_order(&11.into(), &32942478.into()).unwrap(), 2745206);
    }

    #[test]
    fn primitive_roots() {
        // The order of a primitive root is the totient of the modulus.
        for (g, n, totient) in [
            ("7", "10", "4"),
            ("2", "9", "6"),
            ("47", "82", "40"),
            ("3", "686", "294"),
            ("5", "97", "96"),
            ("5", "9409", "9312"),
            ("5", "40487", "40486"),
            ("5", "80974", "40486"),
            ("40492", "1639197169", "1639156682"),
            ("2", "12157665459056928801", "8105110306037952534"),
        ] {
            assert_eq!(n_order(&int(g), &int(n)).unwrap(), int(totient), "n = {n}");
        }
    }

    #[test]
    fn large_prime_power() {
        let p = Integer::from(10).pow(50u32) + 151u32;
        let n = p.clone().square();
        let expected = int(
            "10000000000000000000000000000000000000000000000030100000000000000000000000000000000000000000000022650",
        );
        assert_eq!(n_order(&11.into(), &n).unwrap(), expected);
        assert_eq!(
            n_order_with_factors(&11.into(), &n, &HashMap::from([(p, 2)])).unwrap(),
            expected
        );
    }

    #[test]
    fn several_large_prime_factors() {
        // Both prime factors are out of reach of trial division: the order can only be right if
        // the modulus is really factored.
        let n = Integer::from(33554467u32) * 33554473u32;
        assert_eq!(n_order(&3.into(), &n).unwrap(), int("5212511064222"));
    }

    #[test]
    fn trivial_case() {
        // a % n == 1
        assert_eq!(n_order(&1.into(), &7.into()).unwrap(), 1);
        assert_eq!(n_order(&8.into(), &7.into()).unwrap(), 1);
        assert_eq!(n_order(&15.into(), &7.into()).unwrap(), 1);
        // Negative `a` is reduced modulo `n` first.
        assert_eq!(n_order(&(-6).into(), &7.into()).unwrap(), 1);
        assert_eq!(n_order(&(-1).into(), &7.into()).unwrap(), 2);
    }

    #[test]
    fn not_relatively_prime() {
        assert_eq!(
            n_order(&6.into(), &9.into()),
            Err(Error::NotRelativelyPrime)
        );
        assert_eq!(
            n_order(&0.into(), &7.into()),
            Err(Error::NotRelativelyPrime)
        );
    }

    #[test]
    fn lambda_or_totient_as_the_bound() {
        // `lambda(32942478) = 5490412` where `totient(32942478) = 10980824`: the order of a unit is
        // stripped from the smaller of the two, and a base sharing a factor with the modulus keeps
        // the larger, the only one that bounds its powers.
        let n = int("32942478");
        let n_factors = fast_factor(&n);

        let (order, order_factors) = element_order_with_factors(&n, &11.into(), &n_factors);
        assert_eq!(order, 2745206);
        assert_eq!(order, n_order(&11.into(), &n).unwrap());
        assert_eq!(order_factors, fast_factor(&order));

        // 6 shares the factors 2 and 3 with `n`: no power of it is 1, nothing is stripped, and the
        // totient comes back as it is.
        let (order, order_factors) = element_order_with_factors(&n, &6.into(), &n_factors);
        assert_eq!(order, 10980824);
        assert_eq!(order_factors, fast_factor(&order));
    }

    #[test]
    fn element_order_over_every_small_modulus() {
        // Every base modulo every small modulus: a unit gets its real order, and a base sharing a
        // factor with the modulus gets the totient, which lambda must not replace. The logarithm
        // of 4 in base 2 modulo 12 is 2, and the callers search `0..order` for it: `lambda(12)` is
        // 2 and would stop before reaching it, where `totient(12)` is 4.
        for n in 2u32..100 {
            let n = Integer::from(n);
            let n_factors = fast_factor(&n);
            for b in 0..n.to_u32().unwrap() {
                let b = Integer::from(b);
                let (order, order_factors) = element_order_with_factors(&n, &b, &n_factors);
                // The factorization returned is always the factorization of the order.
                assert_eq!(
                    order_factors,
                    fast_factor(&order),
                    "order {order} of {b} modulo {n}"
                );

                if Integer::from(b.gcd_ref(&n)) == 1 {
                    assert_eq!(order, n_order(&b, &n).unwrap(), "{b} modulo {n}");
                } else {
                    assert_eq!(order, totient(&n_factors), "{b} modulo {n}");
                }
            }
        }
        assert_eq!(
            element_order_with_factors(&12.into(), &2.into(), &fast_factor(&12.into())).0,
            4
        );
    }

    #[test]
    fn exponent_of_the_units_modulo_a_power_of_two() {
        assert_eq!(two_power_exponent(1), 0);
        assert_eq!(two_power_exponent(2), 1);
        for e in 3..64 {
            assert_eq!(two_power_exponent(e), e - 2);
        }

        for e in 1..12u32 {
            let n = Integer::from(1) << e;
            let lambda = Integer::from(1) << two_power_exponent(e as usize) as u32;
            // It really is an exponent of the group: every odd residue is 1 to it.
            for b in (1..n.to_u32().unwrap()).step_by(2) {
                assert_eq!(
                    Integer::from(b).pow_mod(&lambda, &n).unwrap(),
                    1,
                    "{b}**{lambda} mod {n}"
                );
            }
            // And the smallest one, from `e = 3` on: 3 reaches it.
            if e >= 3 {
                assert_eq!(n_order(&3.into(), &n).unwrap(), lambda, "n = {n}");
            }
        }
    }

    #[test]
    fn powers_of_two() {
        // A unit modulo `2**e` gets its order, and a base that is not one keeps the totient: the
        // exponent of the units is half of it from `e = 3` on, and only the unit may start there.
        for e in 1..20u32 {
            let n = Integer::from(1) << e;
            let n_factors = fast_factor(&n);
            for b in 0..16.min(n.to_u32().unwrap()) {
                let b = Integer::from(b);
                let (order, order_factors) = element_order_with_factors(&n, &b, &n_factors);
                assert_eq!(
                    order_factors,
                    fast_factor(&order),
                    "order {order} of {b} modulo {n}"
                );
                if Integer::from(b.gcd_ref(&n)) == 1 {
                    assert_eq!(order, n_order(&b, &n).unwrap(), "{b} modulo {n}");
                } else {
                    assert_eq!(order, totient(&n_factors), "{b} modulo {n}");
                }
            }
        }

        // A power of two mixed with an odd prime, where lambda is an lcm: `lambda(24)` is
        // `lcm(2, 2) = 2`, where the order of the units modulo 8 alone is 4.
        let n = Integer::from(24);
        assert_eq!(
            element_order_with_factors(&n, &5.into(), &fast_factor(&n)).0,
            2
        );
        // The totient is still what a base sharing a factor with the modulus gets: the logarithm of
        // 4 in base 2 modulo 12 is 2, and the callers search `0..order` for it.
        let n = Integer::from(12);
        assert_eq!(
            element_order_with_factors(&n, &2.into(), &fast_factor(&n)).0,
            4
        );
    }

    /// Euler's totient of the number `n_factors` is the factorization of.
    fn totient(n_factors: &HashMap<Integer, usize>) -> Integer {
        n_factors.iter().fold(Integer::from(1), |totient, (p, e)| {
            totient * Integer::from(p - 1u32) * p.clone().pow(*e as u32 - 1)
        })
    }

    #[test]
    fn factorization_validation() {
        // Every one of these used to panic on a subtraction that overflowed, a negative exponent, a
        // division by zero, or run for minutes, or never end at all.
        for (a, n, factors) in [
            // An exponent of 0, which `kx - 1` and `pow(e - 1)` both read as a huge one.
            (2, 9, vec![(3, 0)]),
            (2, 9, vec![(3, usize::MAX)]),
            // Primes nothing is a unit modulo.
            (2, 9, vec![(0, 2)]),
            (3, 4, vec![(0, 3)]),
            (2, 9, vec![(1, 1)]),
            (2, 9, vec![(1, 2)]),
            // Negative primes.
            (3, 4, vec![(-3, 2)]),
            (3, 4, vec![(-5, 2)]),
            // The factorization of another number, whose prime powers multiply to something else:
            // the loop of `prime_power_order` never reached 1 and grew the order for ever.
            (3, 4, vec![(3, 1)]),
            (7, 9, vec![(7, 1)]),
            (2, 9, vec![(3, 1), (5, 2)]),
            (2, 9, vec![(3, 1)]),
            (2, 9, vec![(3, 3)]),
            // A key that is not prime, which the product check cannot see: the bounded loop of
            // `prime_power_order` is what refuses it.
            (3, 4, vec![(4, 1)]),
            (2, 9, vec![(9, 1)]),
        ] {
            let n_factors: HashMap<Integer, usize> = factors
                .iter()
                .map(|&(p, e)| (Integer::from(p), e))
                .collect();
            assert_eq!(
                n_order_with_factors(&a.into(), &n.into(), &n_factors),
                Err(Error::InvalidOrder),
                "a = {a}, n = {n}, factors = {factors:?}"
            );
        }

        // And the factorization that does belong to `n` still gives the order.
        assert_eq!(
            n_order_with_factors(
                &2.into(),
                &9.into(),
                &HashMap::from([(Integer::from(3), 2)])
            )
            .unwrap(),
            6
        );
        // A base that is not relatively prime with `n` is still refused as such, the factorization
        // being checked first.
        assert_eq!(
            n_order_with_factors(
                &3.into(),
                &9.into(),
                &HashMap::from([(Integer::from(3), 2)])
            ),
            Err(Error::NotRelativelyPrime)
        );
    }

    #[test]
    fn modulus_validation() {
        // Modulo 1, every number is 1.
        assert_eq!(n_order(&2.into(), &1.into()).unwrap(), 1);
        assert_eq!(n_order(&0.into(), &1.into()).unwrap(), 1);

        assert_eq!(n_order(&2.into(), &0.into()), Err(Error::InvalidModulus));
        assert_eq!(n_order(&2.into(), &(-1).into()), Err(Error::InvalidModulus));
        assert_eq!(
            n_order_with_factors(&2.into(), &0.into(), &HashMap::new()),
            Err(Error::InvalidModulus)
        );
    }
}
