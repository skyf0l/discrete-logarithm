use rug::Integer;

use crate::Error;

/// Trial multiplication algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// The algorithm finds the discrete logarithm using exhaustive search.
/// This naive method is used as fallback algorithm of ``discrete_log`` when the group order is very small.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_trial_mul(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);
    let order = order.unwrap_or(n);

    // An order that does not fit in a `usize` is out of reach anyway.
    let steps = order.to_usize().unwrap_or(usize::MAX);
    let mut x = Integer::from(1);
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
}
