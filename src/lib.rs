#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use std::collections::HashMap;

use n_order::n_order_with_factors;
use rug::{integer::IsPrime, Integer};
mod n_order;
mod pohlig_hellman;
mod pollard_rho;
mod shanks_steps;
mod trial_mul;
mod utils;

pub use n_order::n_order;
pub use pohlig_hellman::discrete_log_pohlig_hellman;
pub use pollard_rho::discrete_log_pollard_rho;
pub use shanks_steps::discrete_log_shanks_steps;
pub use trial_mul::discrete_log_trial_mul;

/// Discrete logarithm error
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// Log does not exist
    #[error("Log does not exist")]
    LogDoesNotExist,
    /// A and n are not relatively prime
    #[error("A and n are not relatively prime")]
    NotRelativelyPrime,
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
pub fn discrete_log(n: &Integer, a: &Integer, b: &Integer) -> Result<Integer, Error> {
    discrete_log_with_order(n, a, b, &n_order(b, n)?)
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// If the prime factorization of `n` is known, it can be passed as `n_factors` to speed up the computation.
pub fn discrete_log_with_factors(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    discrete_log_with_order(n, a, b, &n_order_with_factors(b, n, n_factors)?)
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_with_order(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
) -> Result<Integer, Error> {
    if *order < 1000 {
        discrete_log_trial_mul(n, a, b, Some(order))
    } else if order.is_probably_prime(100) != IsPrime::No {
        if *order < shanks_steps::MAX_ORDER {
            discrete_log_shanks_steps(n, a, b, Some(order))
        } else {
            discrete_log_pollard_rho(n, a, b, Some(order))
        }
    } else {
        discrete_log_pohlig_hellman(n, a, b, Some(order))
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn big_discrete_log() {
        let n = Integer::from_str("83408372012221120677052349409462320990177094246143674474872152829440524098582262384066400107950985845255268335597502228206679771838750219696329523257176739436871327238322817403970284015587320158034304282786944710043150568360761457471641695390427267786485448748458445872307883254297662715749746270343116946519").unwrap();
        let a = Integer::from_str("109770827223661560471527567179288748906402603483328748683689436879660543465776899146036833470531024202351087008847594392666852763100570391337823820240726499421306887565697452868723849092658743267256316770223643723095601213088336064635680075206929620159782416078143076506249031972043819429093074684182845530529249907297736582589125917235222921623698038868900282049587768700860009877737045693722732170123306528145661683416808514556360429554775212088169626620488741903267154641722293484797745665402402381445609873333905772582972140944493849645600529147490903067975300304532955461710562911203871840101407995813072692212").unwrap();
        let b = Integer::from_str("65537").unwrap();

        assert_eq!(
            discrete_log(&n, &a, &b).unwrap(),
            Integer::from_str("495604594360692646132957963901411709").unwrap(),
        );
    }
}
