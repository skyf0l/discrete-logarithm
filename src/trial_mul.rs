use rug::Integer;

use crate::{Error, check_modulus};

/// Trial multiplication algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// The algorithm finds the discrete logarithm using exhaustive search.
/// This naive method is used as fallback algorithm of ``discrete_log`` when the group order is very small.
///
/// A modulus that is not positive is refused with [`Error::InvalidModulus`]. Modulo 1 every residue
/// is 0, so the logarithm is 0.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
/// Without it the whole modulus bounds the search, which is only a search worth waiting for on a
/// small modulus: an order at or above `2**64`, and a modulus at or above it when no order is
/// given, are beyond what this can walk through.
pub fn discrete_log_trial_mul(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    // Modulo 1 every residue is 0, so every exponent is a logarithm and 0 is the smallest one,
    // whatever `order` says: the search below would answer the same for any order of at least 1.
    if *n == 1 {
        return Ok(Integer::new());
    }
    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);
    let order = order.unwrap_or(n);

    // An order that does not fit in a `usize` is out of reach anyway, and nothing is below a
    // negative one: an order given as one is an order with no exponent at all.
    let steps = if *order < 0 {
        0
    } else {
        order.to_usize().unwrap_or(usize::MAX)
    };
    // Reduced modulo `n`, which only matters modulo 1, where 1 is 0 as everything else is.
    let mut x = Integer::from(1) % n;
    for i in 0..steps {
        if x == a {
            return Ok(Integer::from(i));
        }
        x = x * &b % n;
    }

    Err(Error::LogDoesNotExist)
}

#[cfg(test)]
mod tests {
    use rug::ops::Pow;

    use super::*;

    #[test]
    fn trial_mul() {
        assert_eq!(
            discrete_log_trial_mul(&587.into(), &(Integer::from(2).pow(7)), &2.into(), None)
                .unwrap(),
            7
        );
        assert_eq!(
            discrete_log_trial_mul(&941.into(), &(Integer::from(7).pow(18)), &7.into(), None)
                .unwrap(),
            18
        );
        assert_eq!(
            discrete_log_trial_mul(&389.into(), &(Integer::from(3).pow(81)), &3.into(), None)
                .unwrap(),
            81
        );
        assert_eq!(
            discrete_log_trial_mul(&191.into(), &(Integer::from(19).pow(123)), &19.into(), None)
                .unwrap(),
            123
        );
        // The order of the base bounds the search.
        assert_eq!(
            discrete_log_trial_mul(&587.into(), &128.into(), &2.into(), Some(&293.into())).unwrap(),
            7
        );
        assert_eq!(
            discrete_log_trial_mul(&587.into(), &128.into(), &2.into(), Some(&7.into())),
            Err(Error::LogDoesNotExist)
        );
        // Bases that are not relatively prime with the modulus are searched all the same.
        assert_eq!(
            discrete_log_trial_mul(&9.into(), &0.into(), &3.into(), None).unwrap(),
            2
        );
        assert_eq!(
            discrete_log_trial_mul(&10.into(), &3.into(), &2.into(), None),
            Err(Error::LogDoesNotExist)
        );
        // Negative operands are reduced modulo `n`.
        assert_eq!(
            discrete_log_trial_mul(&587.into(), &(-459).into(), &(-585).into(), None).unwrap(),
            7
        );
    }

    #[test]
    fn invalid_modulus() {
        // Nothing is defined modulo these, and reducing by them is not either.
        for n in [0, -1, -4] {
            assert_eq!(
                discrete_log_trial_mul(&n.into(), &1.into(), &3.into(), None),
                Err(Error::InvalidModulus)
            );
            assert_eq!(
                discrete_log_trial_mul(&n.into(), &1.into(), &3.into(), Some(&10.into())),
                Err(Error::InvalidModulus)
            );
        }
    }

    #[test]
    fn modulus_of_one() {
        // Modulo 1 every residue is 0, so the logarithm is 0 whatever `a` and `b` are.
        assert_eq!(
            discrete_log_trial_mul(&1.into(), &0.into(), &0.into(), None).unwrap(),
            0
        );
        assert_eq!(
            discrete_log_trial_mul(&1.into(), &5.into(), &7.into(), Some(&1.into())).unwrap(),
            0
        );
        // Even an order no exponent is below: modulo 1 there is nothing to search for.
        assert_eq!(
            discrete_log_trial_mul(&1.into(), &5.into(), &7.into(), Some(&(-1).into())).unwrap(),
            0
        );
    }

    #[test]
    fn negative_order() {
        // No exponent is below a negative order, `0` included: the search is empty instead of
        // walking the whole of a `usize`.
        for order in [-1, -5] {
            assert_eq!(
                discrete_log_trial_mul(&587.into(), &1.into(), &2.into(), Some(&order.into())),
                Err(Error::LogDoesNotExist)
            );
        }
    }
}
