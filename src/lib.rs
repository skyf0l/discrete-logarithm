#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use std::collections::HashMap;

use n_order::n_order_with_factors;
use rug::Integer;
mod n_order;
mod trial_mul;

pub use n_order::n_order;
pub use trial_mul::{discrete_log_trial_mul, discrete_log_trial_mul_with_order};

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
    discrete_log_trial_mul_with_order(n, a, b, &n_order(b, n)?)
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
    discrete_log_trial_mul_with_order(n, a, b, &n_order_with_factors(b, n, n_factors)?)
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
    discrete_log_trial_mul_with_order(n, a, b, order)
}
