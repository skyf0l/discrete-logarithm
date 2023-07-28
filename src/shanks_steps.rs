use rug::Integer;

use crate::Error;

/// Trial multiplication algorithm for computing the discrete logarithm of `h` in base `g` modulo `n` (smallest non-negative integer `x` where `g**x = h (mod n)`)
///
/// The algorithm finds the discrete logarithm using exhaustive search.
/// This naive method is used as fallback algorithm of ``discrete_log`` when the group order is very small.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
pub fn discrete_log_shanks_steps(
    _n: &Integer,
    _a: &Integer,
    _b: &Integer,
    _order: Option<&Integer>,
) -> Result<Integer, Error> {
    unimplemented!()
}
