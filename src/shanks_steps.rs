use std::collections::HashMap;

use rug::Integer;

use crate::{n_order, Error};

pub const MAX_ORDER: u64 = 1_000_000_000_000u64;

/// Baby-step giant-step algorithm for computing the discrete logarithm of `h` in base `g` modulo `n` (smallest non-negative integer `x` where `g**x = h (mod n)`)
///
/// The algorithm is a time-memory trade-off of the method of exhaustive search. It uses `O(sqrt(m))` memory, where `m` is the group order.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_shanks_steps(
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

    if order >= MAX_ORDER {
        return Err(Error::LogDoesNotExist);
    }

    let m = order.sqrt().to_i32().unwrap() + 1;
    let mut t = HashMap::new();
    let mut x = Integer::from(1);
    for i in 0..m {
        t.insert(x.clone(), i);
        x = x * &b % n;
    }

    let z = b.invert(n).unwrap();
    let z = z.pow_mod(&Integer::from(m), n).unwrap();
    let mut x = a;
    for i in 0..m {
        if let Some(j) = t.get(&x) {
            return Ok(Integer::from(i * m + j));
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
}
