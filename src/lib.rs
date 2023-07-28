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
