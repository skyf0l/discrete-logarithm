#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use rug::Integer;
mod trial_mul;

pub use trial_mul::discrete_log_trial_mul;
use trial_mul::discrete_log_trial_mul_with_order;

/// Discrete logarithm error
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// Log does not exist
    #[error("Log does not exist")]
    LogDoesNotExist,
}

/// Compute the discrete logarithm of `h` in base `g` modulo `n` (smallest non-negative integer `x` where `g**x = h (mod n)`).
pub fn discrete_log(n: &Integer, a: &Integer, b: &Integer) -> Result<Integer, Error> {
    discrete_log_trial_mul_with_order(n, a, b, &None)
}

/// Compute the discrete logarithm of `h` in base `g` modulo `n` (smallest non-negative integer `x` where `g**x = h (mod n)`).
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_with_order(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Option<Integer>,
) -> Result<Integer, Error> {
    discrete_log_trial_mul_with_order(n, a, b, order)
}
