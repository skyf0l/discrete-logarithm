use std::collections::HashMap;

use primal::Primes;
use rug::{ops::Pow, Integer};

use crate::Error;

fn fast_factor(n: &Integer) -> HashMap<Integer, usize> {
    let mut factors: HashMap<Integer, usize> = HashMap::new();
    let mut n: Integer = n.clone();
    for prime in Primes::all().take(1_000_000) {
        let prime = Integer::from(prime);
        if n.clone().div_rem(prime.clone()).1 == 0 {
            // factors.insert(prime.clone(), 1);
            while n.clone().div_rem(prime.clone()).1 == 0 {
                n /= &prime;
                *factors.entry(prime.clone()).or_insert(0) += 1;
            }
        }
    }

    if n != 1 {
        *factors.entry(n).or_insert(0) += 1;
    }

    factors
}

/// Returns the order of `a` modulo `n`.
///
/// The order of `a` modulo `n` is the smallest integer `k` such that `a**k` leaves a remainder of 1 with `n`.
pub fn n_order(a: &Integer, n: &Integer) -> Result<Integer, Error> {
    if a.clone().gcd(n) != 1 {
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
    if a.clone().gcd(n) != 1 {
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
    let a = a.clone() % n;
    for (p, e) in factors {
        let mut exponent = group_order.clone();
        for f in 0..=e {
            if a.clone().pow_mod(&exponent, n).unwrap() != 1 {
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

        assert_eq!(n_order_with_factors(&11.into(), &Integer::from(Integer::from(10).pow(50) + 151).square(), &HashMap::from([(Integer::from(Integer::from(10).pow(50) + 151), 2)])).unwrap(), Integer::from_str("10000000000000000000000000000000000000000000000030100000000000000000000000000000000000000000000022650").unwrap());
    }
}
