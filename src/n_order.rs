use std::collections::HashMap;

use rug::{ops::Pow, Integer};

use crate::{check_modulus, factor::fast_factor, Error};

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
    n_order_with_factors(a, n, &fast_factor(n))
}

/// Returns the order of `a` modulo `n`.
///
/// The order of `a` modulo `n` is the smallest integer `k` such that `a**k` leaves a remainder of 1 with `n`.
///
/// If the prime factorization of `n` is known, it can be passed as `n_factors` to speed up the computation.
pub fn n_order_with_factors(
    a: &Integer,
    n: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::from(1));
    }

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
        order.lcm_mut(&prime_power_order(&a, p, *e));
    }
    Ok(order)
}

/// Returns the order of `a` modulo `p**e`, `a` and `p` being relatively prime.
fn prime_power_order(a: &Integer, p: &Integer, e: usize) -> Integer {
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
        while x != 1 {
            x = x.pow_mod(px, &pe).unwrap();
            order *= px;
        }
    }
    order
}

/// Returns the order of `b` modulo `n` and the prime factorization of that order.
///
/// The order of the whole group of units and its factorization are computed in one pass, then
/// the primes the order of `b` does not need are divided out: the factorization of the order
/// comes for free, and the primality test of [`n_order`] is not needed.
///
/// `b` must be reduced modulo `n`. When `b` and `n` are not relatively prime, the powers of `b`
/// are not a subgroup of the units and nothing can be divided out: the returned order is the
/// order of the group, which the callers only use as a bound.
pub fn element_order_with_factors(
    n: &Integer,
    b: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> (Integer, HashMap<Integer, usize>) {
    // The order of the group of units is `totient(n)`, computed together with its factorization.
    let mut factors: HashMap<Integer, usize> = HashMap::new();
    for (px, kx) in n_factors {
        if *kx > 1 {
            *factors.entry(px.clone()).or_insert(0) += kx - 1;
        }
        for (py, ky) in fast_factor(&Integer::from(px - 1u32)) {
            *factors.entry(py).or_insert(0) += ky;
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
        let expected = int("10000000000000000000000000000000000000000000000030100000000000000000000000000000000000000000000022650");
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
