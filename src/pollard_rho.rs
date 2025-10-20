use crate::bignum::{Integer, new_rng};

use crate::{n_order, Error};

const RETRIES: usize = 10;

/// Pollard's Rho  algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// It is a randomized algorithm with the same expected running time as `discrete_log_shanks_steps`, but requires a negligible amount of memory.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_pollard_rho(
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

    let mut rand_state = new_rng();

    let order_minus_2 = Integer::from(&order - 2);

    for _ in 0..RETRIES {
        let mut aa = order_minus_2.clone().random_below(&mut rand_state) + 1;
        let mut ba = order_minus_2.clone().random_below(&mut rand_state) + 1;
        let mut xa = b.clone().pow_mod(&aa, n).unwrap() * a.clone().pow_mod(&ba, n).unwrap() % n;

        let c = xa.clone() % 3;
        let mut xb;
        let mut ab;
        let mut bb;
        if c == 0 {
            xb = a.clone() * &xa % n;
            ab = aa.clone();
            bb = (ba.clone() + 1) % &order;
        } else if c == 1 {
            xb = xa.clone() * &xa % n;
            ab = (aa.clone() + &aa) % &order;
            bb = (ba.clone() + &ba) % &order;
        } else {
            xb = b.clone() * &xa % n;
            ab = (aa.clone() + 1) % &order;
            bb = ba.clone();
        }

        for _ in 0..order.to_u32().unwrap_or(u32::MAX) {
            let c = xa.clone() % 3;
            if c == 0 {
                xa = a.clone() * &xa % n;
                ba = (ba.clone() + 1) % &order;
            } else if c == 1 {
                xa = xa.clone() * &xa % n;
                aa = (aa.clone() + &aa) % &order;
                ba = (ba.clone() + &ba) % &order;
            } else {
                xa = b.clone() * &xa % n;
                aa = (aa.clone() + 1) % &order;
            }

            let c = xb.clone() % 3;
            if c == 0 {
                xb = a.clone() * &xb % n;
                bb = (bb.clone() + 1) % &order;
            } else if c == 1 {
                xb = xb.clone() * &xb % n;
                ab = (ab.clone() + &ab) % &order;
                bb = (bb.clone() + &bb) % &order;
            } else {
                xb = b.clone() * &xb % n;
                ab = (ab.clone() + 1) % &order;
            }

            let c = xb.clone() % 3;
            if c == 0 {
                xb = a.clone() * &xb % n;
                bb = (bb.clone() + 1) % &order;
            } else if c == 1 {
                xb = xb.clone() * &xb % n;
                ab = (ab.clone() + &ab) % &order;
                bb = (bb.clone() + &bb) % &order;
            } else {
                xb = b.clone() * &xb % n;
                ab = (ab.clone() + 1) % &order;
            }

            if xa == xb {
                let r = (ba.clone() - &bb) % &order;
                if let Ok(i) = r.invert(&order) {
                    let e = (i * (ab.clone() - aa.clone()) % &order + &order) % &order;
                    if (b.clone().pow_mod(&e, n).unwrap() - &a) % n == 0 {
                        return Ok(e);
                    }
                }
                break;
            }
        }
    }

    Err(Error::LogDoesNotExist)
}

#[cfg(test)]
mod tests {
    use rug::ops::Pow;

    use super::*;

    #[test]
    fn pollard_rho() {
        assert_eq!(
            discrete_log_pollard_rho(&6013199.into(), &(Integer::from(2).pow(6)), &2.into(), None)
                .unwrap(),
            6
        );
        assert_eq!(
            discrete_log_pollard_rho(
                &6138719.into(),
                &(Integer::from(2).pow(19)),
                &2.into(),
                None
            )
            .unwrap(),
            19
        );
        assert_eq!(
            discrete_log_pollard_rho(
                &36721943.into(),
                &(Integer::from(2).pow(40)),
                &2.into(),
                None
            )
            .unwrap(),
            40
        );
        assert_eq!(
            discrete_log_pollard_rho(
                &24567899.into(),
                &(Integer::from(3).pow(333)),
                &3.into(),
                None
            )
            .unwrap(),
            333
        );
        assert_eq!(
            discrete_log_pollard_rho(&11.into(), &7.into(), &31.into(), None),
            Err(Error::LogDoesNotExist)
        );
        assert_eq!(
            discrete_log_pollard_rho(&227.into(), &(Integer::from(3).pow(7)), &5.into(), None),
            Err(Error::LogDoesNotExist)
        );
    }
}
