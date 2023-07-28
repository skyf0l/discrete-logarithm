#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use rug::Integer;

/// Compute the discrete logarithm of `h` in base `g` modulo `n`.
/// Returns the smallest non-negative integer `x` where `g**x = h (mod n)`.
pub fn discrete_log(_n: &Integer, _a: &Integer, _b: &Integer) -> Integer {
    todo!()
}
