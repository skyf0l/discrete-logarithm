use std::collections::HashMap;

use rug::{ops::Pow, Integer};

use crate::{utils::fast_factor, Error};

/// Returns the order of `a` modulo `n`.
///
/// The order of `a` modulo `n` is the smallest integer `k` such that `a**k` leaves a remainder of 1 with `n`.
pub fn n_order(a: &Integer, n: &Integer) -> Result<Integer, Error> {
    // Validate n > 1
    if *n <= 1 {
        return Err(Error::NotRelativelyPrime);
    }

    // Early return for trivial case
    let a_mod = a.clone() % n;
    if a_mod == 1 {
        return Ok(Integer::from(1));
    }

    if a_mod.clone().gcd(n) != 1 {
        return Err(Error::NotRelativelyPrime);
    }

    let factors = fast_factor(n);
    n_order_with_factors(a, n, &factors)
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
    // Validate n > 1
    if *n <= 1 {
        return Err(Error::NotRelativelyPrime);
    }

    // Early return for trivial case
    let a_mod = a.clone() % n;
    if a_mod == 1 {
        return Ok(Integer::from(1));
    }

    if a_mod.clone().gcd(n) != 1 {
        return Err(Error::NotRelativelyPrime);
    }

    let mut factors = HashMap::new();
    for (px, kx) in n_factors.iter() {
        if *kx > 1 {
            *factors.entry(px.clone()).or_insert(0) += kx - 1;
        }
        let fpx = fast_factor(&(px.clone() - 1));
        for (py, ky) in fpx.iter() {
            *factors.entry(py.clone()).or_insert(0) += ky;
        }
    }

    let mut group_order = Integer::from(1);
    for (px, kx) in factors.iter() {
        group_order *= px.clone().pow(*kx as u32);
    }

    let mut order = Integer::from(1);
    for (p, e) in factors {
        let mut exponent = group_order.clone();
        for f in 0..=e {
            if a_mod.clone().pow_mod(&exponent, n).unwrap() != 1 {
                order *= p.clone().pow((e - f + 1) as u32);
                break;
            }
            exponent /= &p;
        }
    }

    Ok(order)
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn trial_mul() {
        assert_eq!(n_order(&2.into(), &13.into()).unwrap(), 12);
        for (a, res) in (1..=6).zip(vec![1, 3, 6, 3, 6, 2]) {
            assert_eq!(n_order(&a.into(), &7.into()).unwrap(), res);
        }
        assert_eq!(n_order(&5.into(), &17.into()).unwrap(), 16);
        assert_eq!(
            n_order(&17.into(), &11.into()),
            n_order(&6.into(), &11.into())
        );
        assert_eq!(n_order(&101.into(), &119.into()).unwrap(), 6);
        assert_eq!(
            n_order(&6.into(), &9.into()),
            Err(Error::NotRelativelyPrime)
        );

        assert_eq!(n_order_with_factors(&11.into(), &(Integer::from(10).pow(50) + 151u64).square(), &HashMap::from([(Integer::from(10).pow(50) + 151, 2)])).unwrap(), Integer::from_str("10000000000000000000000000000000000000000000000030100000000000000000000000000000000000000000000022650").unwrap());
    }

    #[test]
    fn n_order_trivial_case() {
        // Test early return for a % n == 1
        assert_eq!(n_order(&1.into(), &7.into()).unwrap(), 1);
        assert_eq!(n_order(&8.into(), &7.into()).unwrap(), 1); // 8 % 7 = 1
        assert_eq!(n_order(&15.into(), &7.into()).unwrap(), 1); // 15 % 7 = 1
    }

    #[test]
    fn n_order_validation() {
        // Test n <= 1 validation
        assert_eq!(
            n_order(&2.into(), &1.into()),
            Err(Error::NotRelativelyPrime)
        );
        assert_eq!(
            n_order(&2.into(), &0.into()),
            Err(Error::NotRelativelyPrime)
        );
        assert_eq!(
            n_order(&2.into(), &(-1).into()),
            Err(Error::NotRelativelyPrime)
        );
    }
}
