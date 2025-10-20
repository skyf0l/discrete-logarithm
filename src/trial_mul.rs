use crate::bignum::Integer;

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
    let a = a.clone() % n;
    let b = b.clone() % n;
    let order = match order {
        Some(order) => order,
        None => n,
    };

    let mut x = Integer::from(1);
    let mut i = 0;
    loop {
        if x == a {
            return Ok(Integer::from(i));
        }
        x = x * &b % n;

        i += 1;
        if i == *order {
            break;
        }
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
    }
}
