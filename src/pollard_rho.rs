use rug::{rand::RandState, Integer};

use crate::{discrete_log_trial_mul, n_order, Error};

/// Number of random starting points tried before giving up.
const RETRIES: usize = 10;

/// Pollard's Rho algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
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
    solve(n, a, b, order, &mut RandState::new())
}

/// Pollard's Rho algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// Same as [`discrete_log_pollard_rho`], with the random generator seeded with `seed`: the same
/// seed always tries the same starting points, and a retry with another seed makes other choices.
pub fn discrete_log_pollard_rho_with_seed(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    seed: u64,
) -> Result<Integer, Error> {
    let mut rand_state = RandState::new();
    rand_state.seed(&Integer::from(seed));
    solve(n, a, b, order, &mut rand_state)
}

fn solve(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error> {
    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);
    let order = match order {
        Some(order) => order.clone(),
        None => n_order(&b, n)?,
    };

    // There is no room for random starting points, and an exhaustive search costs nothing.
    if order < 4 {
        return discrete_log_trial_mul(n, &a, &b, Some(&order));
    }

    // Exponents are drawn in `[1, order - 1]`.
    let order_minus_1 = Integer::from(&order - 1u32);
    // An order that does not fit in a `usize` is out of reach anyway.
    let steps = order.to_usize().unwrap_or(usize::MAX);

    for _ in 0..RETRIES {
        let mut aa = order_minus_1.clone().random_below(rand_state) + 1u32;
        let mut ba = order_minus_1.clone().random_below(rand_state) + 1u32;
        let mut xa = b.clone().pow_mod(&aa, n).unwrap() * a.clone().pow_mod(&ba, n).unwrap() % n;

        let c = xa.clone() % 3;
        let mut xb;
        let mut ab;
        let mut bb;
        if c == 0 {
            xb = a.clone() * &xa % n;
            ab = aa.clone();
            bb = (ba.clone() + 1u32) % &order;
        } else if c == 1 {
            xb = xa.clone() * &xa % n;
            ab = (aa.clone() + &aa) % &order;
            bb = (ba.clone() + &ba) % &order;
        } else {
            xb = b.clone() * &xa % n;
            ab = (aa.clone() + 1u32) % &order;
            bb = ba.clone();
        }

        for _ in 0..steps {
            let c = xa.clone() % 3;
            if c == 0 {
                xa = a.clone() * &xa % n;
                ba = (ba.clone() + 1u32) % &order;
            } else if c == 1 {
                xa = xa.clone() * &xa % n;
                aa = (aa.clone() + &aa) % &order;
                ba = (ba.clone() + &ba) % &order;
            } else {
                xa = b.clone() * &xa % n;
                aa = (aa.clone() + 1u32) % &order;
            }

            // The second walker moves twice as fast: the cycle is found without any memory.
            for _ in 0..2 {
                let c = xb.clone() % 3;
                if c == 0 {
                    xb = a.clone() * &xb % n;
                    bb = (bb.clone() + 1u32) % &order;
                } else if c == 1 {
                    xb = xb.clone() * &xb % n;
                    ab = (ab.clone() + &ab) % &order;
                    bb = (bb.clone() + &bb) % &order;
                } else {
                    xb = b.clone() * &xb % n;
                    ab = (ab.clone() + 1u32) % &order;
                }
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

    #[test]
    fn seeded() {
        // The same seed always gives the same walk, and the result does not depend on it.
        for seed in [0, 1, 42] {
            assert_eq!(
                discrete_log_pollard_rho_with_seed(
                    &6013199.into(),
                    &(Integer::from(2).pow(6)),
                    &2.into(),
                    None,
                    seed
                )
                .unwrap(),
                6
            );
            assert_eq!(
                discrete_log_pollard_rho_with_seed(
                    &24567899.into(),
                    &(Integer::from(3).pow(333)),
                    &3.into(),
                    None,
                    seed
                )
                .unwrap(),
                333
            );
            assert_eq!(
                discrete_log_pollard_rho_with_seed(&11.into(), &7.into(), &31.into(), None, seed),
                Err(Error::LogDoesNotExist)
            );
        }
    }

    #[test]
    fn tiny_orders() {
        // Orders too small to draw a starting point from.
        for order in 1..4u32 {
            assert_eq!(
                discrete_log_pollard_rho(&7.into(), &1.into(), &2.into(), Some(&order.into()))
                    .unwrap(),
                0
            );
        }
        assert_eq!(
            discrete_log_pollard_rho(&7.into(), &2.into(), &2.into(), Some(&1.into())),
            Err(Error::LogDoesNotExist)
        );
    }
}
