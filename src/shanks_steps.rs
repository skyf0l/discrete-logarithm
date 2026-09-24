use std::collections::HashMap;

use rug::Integer;

use crate::{n_order, Error};

/// Orders of this size or more need more memory than the table can hold.
pub const MAX_ORDER: u64 = 1_000_000_000_000u64;

/// Baby-step giant-step algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// The algorithm is a time-memory trade-off of the method of exhaustive search. It uses `O(sqrt(m))` memory, where `m` is the group order.
///
/// Orders of 10^12 or more are refused with [`Error::LogDoesNotExist`]: the table of the baby
/// steps would not fit in memory. Use [`discrete_log_pollard_rho`](crate::discrete_log_pollard_rho) for those.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_shanks_steps(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);
    let order = match order {
        Some(order) => order.clone(),
        None => n_order(&b, n)?,
    };

    if order >= MAX_ORDER {
        return Err(Error::LogDoesNotExist);
    }

    let m = order.sqrt() + 1u32;
    // `m` is at most the square root of `MAX_ORDER`: the steps fit in a `u64`.
    let steps = m.to_u64().unwrap();

    // Baby steps: b**i for every i < m.
    let mut table = HashMap::with_capacity(steps as usize);
    let mut x = Integer::from(1);
    for i in 0..steps {
        table.insert(x.clone(), i);
        x = x * &b % n;
    }

    // Giant steps: a * b**(-i*m) for every i < m.
    let z = b
        .invert(n)
        .map_err(|_| Error::NotRelativelyPrime)?
        .pow_mod(&m, n)
        .unwrap();
    let mut x = a;
    for i in 0..steps {
        if let Some(j) = table.get(&x) {
            return Ok(Integer::from(i) * &m + j);
        }
        x = x * &z % n;
    }

    Err(Error::LogDoesNotExist)
}

#[cfg(test)]
mod tests {
    use rug::ops::Pow;

    use super::*;

    #[test]
    fn shanks_steps() {
        assert_eq!(
            discrete_log_shanks_steps(&442879.into(), &(Integer::from(7).pow(2)), &7.into(), None)
                .unwrap(),
            2
        );
        assert_eq!(
            discrete_log_shanks_steps(&874323.into(), &(Integer::from(5).pow(19)), &5.into(), None)
                .unwrap(),
            19
        );
        assert_eq!(
            discrete_log_shanks_steps(
                &6876342.into(),
                &(Integer::from(7).pow(71)),
                &7.into(),
                None
            )
            .unwrap(),
            71
        );
        assert_eq!(
            discrete_log_shanks_steps(
                &2456747.into(),
                &(Integer::from(3).pow(321)),
                &3.into(),
                None
            )
            .unwrap(),
            321
        );
    }

    #[test]
    fn no_logarithm() {
        // The order is too large for the table.
        assert_eq!(
            discrete_log_shanks_steps(
                &Integer::from(265390227570863u64),
                &Integer::from(184500076053622u64),
                &2.into(),
                Some(&Integer::from(MAX_ORDER))
            ),
            Err(Error::LogDoesNotExist)
        );
        // The base is not invertible modulo `n`.
        assert_eq!(
            discrete_log_shanks_steps(&32942478.into(), &6.into(), &6.into(), Some(&1000.into())),
            Err(Error::NotRelativelyPrime)
        );
        // No power of the base is `a`.
        assert_eq!(
            discrete_log_shanks_steps(&442879.into(), &3.into(), &7.into(), None),
            Err(Error::LogDoesNotExist)
        );
    }
}
