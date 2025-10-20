use crate::bignum::{Integer, Pow, IntegerExt as _};

use crate::{
    discrete_log_with_order, n_order,
    utils::{crt, fast_factor},
    Error,
};

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
    let a = a.clone() % n;
    let b = b.clone() % n;
    let order = match order {
        Some(order) => order.clone(),
        None => n_order(&b, n)?,
    };

    let order_factors = fast_factor(&order);
    let mut residues = (0..order_factors.len())
        .map(|_| Integer::from(0))
        .collect::<Vec<_>>();

    for (i, (pi, ri)) in order_factors.iter().enumerate() {
        for j in 0..*ri as u32 {
            let gj = b.clone().pow_mod(&residues[i], n).unwrap();
            let aj = (&a * gj.clone().invert(n).unwrap())
                .pow_mod(&(&order / pi.clone().pow(j + 1)), n)
                .unwrap();
            let bj = b.clone().pow_mod(&(&order / pi.clone()), n).unwrap();
            let cj = discrete_log_with_order(n, &aj, &bj, pi)?;
            residues[i] += &cj * pi.clone().pow(j);
        }
    }

    let modulis = order_factors
        .iter()
        .map(|(pi, ri)| pi.clone().pow(*ri as u32))
        .collect::<Vec<_>>();

    if let Some(d) = crt(&residues, &modulis) {
        Ok(d)
    } else {
        Err(Error::LogDoesNotExist)
    }
}

#[cfg(test)]
mod tests {
    use rug::ops::Pow;

    use super::*;

    #[test]
    fn pollard_rho() {
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
}
