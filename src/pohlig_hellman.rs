use std::collections::HashMap;

use rug::{ops::Pow, Integer};

use crate::{crt::crt, discrete_log_with_prime_order, factor::fast_factor, n_order, Error};

/// Pohlig-Hellman algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// In order to compute the discrete logarithm, the algorithm takes advantage of the factorization of the group order. It is more efficient when the group order factors into many small primes.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_pohlig_hellman(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    let order = match order {
        Some(order) => order.clone(),
        None => n_order(b, n)?,
    };
    let order_factors = fast_factor(&order);
    solve(n, a, b, &order, &order_factors)
}

/// Pohlig-Hellman algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// Same as [`discrete_log_pohlig_hellman`] with the prime factorization of the order known: it
/// is the factorization the algorithm is built on, so nothing is left to compute.
pub fn discrete_log_pohlig_hellman_with_factors(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    order_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    solve(n, a, b, order, order_factors)
}

fn solve(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    order_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);

    // One residue per prime power of the order, from a discrete logarithm modulo each prime.
    let mut residues = Vec::with_capacity(order_factors.len());
    let mut modulli = Vec::with_capacity(order_factors.len());

    for (pi, ri) in order_factors {
        let mut residue = Integer::new();
        for j in 0..*ri as u32 {
            let gj = b
                .clone()
                .pow_mod(&residue, n)
                .unwrap()
                .invert(n)
                .map_err(|_| Error::NotRelativelyPrime)?;
            let aj = (&a * gj)
                .pow_mod(&(order.clone() / pi.clone().pow(j + 1)), n)
                .unwrap();
            let bj = b.clone().pow_mod(&(order.clone() / pi.clone()), n).unwrap();
            // `pi` is prime: the sub-problem never comes back to Pohlig-Hellman.
            let cj = discrete_log_with_prime_order(n, &aj, &bj, pi)?;
            residue += cj * pi.clone().pow(j);
        }
        residues.push(residue);
        modulli.push(pi.clone().pow(*ri as u32));
    }

    crt(&residues, &modulli).ok_or(Error::LogDoesNotExist)
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
}
