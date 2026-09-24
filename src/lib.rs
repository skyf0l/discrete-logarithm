#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use std::collections::HashMap;

use rug::{integer::IsPrime, Integer};

use factor::fast_factor;
use n_order::element_order_with_factors;

mod crt;
mod factor;
mod index_calculus;
mod n_order;
mod pohlig_hellman;
mod pollard_rho;
mod shanks_steps;
mod trial_mul;

pub use index_calculus::{discrete_log_index_calculus, discrete_log_index_calculus_with_seed};
pub use n_order::{n_order, n_order_with_factors};
pub use pohlig_hellman::{discrete_log_pohlig_hellman, discrete_log_pohlig_hellman_with_factors};
pub use pollard_rho::{discrete_log_pollard_rho, discrete_log_pollard_rho_with_seed};
pub use shanks_steps::discrete_log_shanks_steps;
pub use trial_mul::discrete_log_trial_mul;

/// Internals exposed for benchmarks only. Not part of the public API, no stability guarantees.
#[cfg(feature = "bench")]
#[doc(hidden)]
pub mod bench {
    pub use crate::factor::fast_factor;
    pub use crate::index_calculus::is_smooth;
    pub use crate::n_order::element_order_with_factors;
}

/// Discrete logarithm error
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// Log does not exist
    #[error("Log does not exist")]
    LogDoesNotExist,
    /// A and n are not relatively prime
    #[error("A and n are not relatively prime")]
    NotRelativelyPrime,
    /// The modulus is not positive
    #[error("n should be positive")]
    InvalidModulus,
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// The order of `b` and its prime factorization are computed in one pass, then the algorithm
/// best suited to that order solves the problem.
pub fn discrete_log(n: &Integer, a: &Integer, b: &Integer) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::new());
    }
    discrete_log_with_factors(n, a, b, &fast_factor(n))
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
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::new());
    }

    let b = b.clone().modulo(n);
    let (order, order_factors) = element_order_with_factors(n, &b, n_factors);
    // A single prime with exponent 1: the order is prime, no primality test needed.
    let prime_order = order_factors.len() == 1 && order_factors.values().all(|&e| e == 1);
    solve(n, a, &b, &order, Some(&order_factors), Some(prime_order))
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
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::new());
    }
    solve(n, a, b, order, None, None)
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// Same as [`discrete_log_with_order`], for an order known to be prime: the primality test of
/// the order is skipped. Passing a composite order gives a wrong result or no result at all.
pub fn discrete_log_with_prime_order(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    if *n == 1 {
        return Ok(Integer::new());
    }
    solve(n, a, b, order, None, Some(true))
}

/// Solves the problem with the algorithm best suited to `order`.
///
/// `order_factors` and `prime_order`, when known, save a factorization and a primality test.
fn solve(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    order_factors: Option<&HashMap<Integer, usize>>,
    prime_order: Option<bool>,
) -> Result<Integer, Error> {
    if *order < 1000 {
        return discrete_log_trial_mul(n, a, b, Some(order));
    }

    let prime_order = prime_order.unwrap_or_else(|| order.is_probably_prime(100) != IsPrime::No);
    if !prime_order {
        return match order_factors {
            Some(order_factors) => {
                discrete_log_pohlig_hellman_with_factors(n, a, b, order, order_factors)
            }
            None => discrete_log_pohlig_hellman(n, a, b, Some(order)),
        };
    }

    // Shanks and Pollard rho are O(sqrt(order)) while index calculus is O(exp(2*sqrt(log(n)log(log(n)))))
    // we compare the expected running times to determine the algorithm which is expected to be faster
    let log_n = n.to_f64().ln();
    let log_log_n = log_n.ln();
    let log_order = order.to_f64().ln();

    // Use index calculus if 4*sqrt(log(n)*log(log(n))) < log(order) - 10
    if 4.0 * (log_n * log_log_n).sqrt() < log_order - 10.0 {
        discrete_log_index_calculus(n, a, b, Some(order))
    } else if *order < shanks_steps::MAX_ORDER {
        // Shanks seems typically faster, but uses O(sqrt(order)) memory
        discrete_log_shanks_steps(n, a, b, Some(order))
    } else {
        discrete_log_pollard_rho(n, a, b, Some(order))
    }
}

/// Rejects the moduli no logarithm and no order are defined for.
pub(crate) fn check_modulus(n: &Integer) -> Result<(), Error> {
    if *n < 1 {
        Err(Error::InvalidModulus)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use rug::ops::Pow;

    use super::*;

    fn int(s: &str) -> Integer {
        Integer::from_str(s).unwrap()
    }

    #[test]
    fn discrete_log_() {
        assert_eq!(
            discrete_log(&587.into(), &(Integer::from(2).pow(9)), &2.into(),).unwrap(),
            9
        );
        assert_eq!(
            discrete_log(&2456747.into(), &(Integer::from(3).pow(51)), &3.into(),).unwrap(),
            51
        );
        assert_eq!(
            discrete_log(&32942478.into(), &(Integer::from(11).pow(127)), &11.into(),).unwrap(),
            127
        );
        assert_eq!(
            discrete_log(
                &int("432751500361"),
                &(Integer::from(7).pow(324)),
                &7.into(),
            )
            .unwrap(),
            324
        );
        assert_eq!(
            discrete_log(&int("265390227570863"), &int("184500076053622"), &2.into(),).unwrap(),
            int("17835221372061"),
        );
        assert_eq!(
            discrete_log(
                &int("22708823198678103974314518195029102158525052496759285596453269189798311427475159776411276642277139650833937"),
                &int("17463946429475485293747680247507700244427944625055089103624311227422110546803452417458985046168310373075327"),
                &123456.into(),
            )
            .unwrap(),
            int("2068031853682195777930683306640554533145512201725884603914601918777510185469769997054750835368413389728895"),
        );
        assert_eq!(
            discrete_log(&5779.into(), &3528.into(), &6215.into(),).unwrap(),
            687
        );
    }

    #[test]
    fn big_discrete_log() {
        let n = int("83408372012221120677052349409462320990177094246143674474872152829440524098582262384066400107950985845255268335597502228206679771838750219696329523257176739436871327238322817403970284015587320158034304282786944710043150568360761457471641695390427267786485448748458445872307883254297662715749746270343116946519");
        let a = int("109770827223661560471527567179288748906402603483328748683689436879660543465776899146036833470531024202351087008847594392666852763100570391337823820240726499421306887565697452868723849092658743267256316770223643723095601213088336064635680075206929620159782416078143076506249031972043819429093074684182845530529249907297736582589125917235222921623698038868900282049587768700860009877737045693722732170123306528145661683416808514556360429554775212088169626620488741903267154641722293484797745665402402381445609873333905772582972140944493849645600529147490903067975300304532955461710562911203871840101407995813072692212");
        let b = int("65537");

        assert_eq!(
            discrete_log(&n, &a, &b).unwrap(),
            int("495604594360692646132957963901411709"),
        );
    }

    #[test]
    fn several_large_prime_factors() {
        // Both prime factors of the modulus are out of reach of trial division: the logarithm can
        // only be found if the modulus is really factored.
        let n = int("1125902456980891");
        assert_eq!(
            discrete_log(&n, &int("55959216458270"), &3.into()).unwrap(),
            123456789
        );
    }

    #[test]
    fn known_order() {
        // The order of 3 modulo 2456747, a prime.
        let n = int("2456747");
        let a = int("406989");
        let order = int("1228373");
        assert_eq!(
            discrete_log_with_order(&n, &a, &3.into(), &order).unwrap(),
            51
        );
        // The same order, known to be prime: one primality test less.
        assert_eq!(
            discrete_log_with_prime_order(&n, &a, &3.into(), &order).unwrap(),
            51
        );
        // A composite order, its factorization unknown: Pohlig-Hellman factors it.
        assert_eq!(
            discrete_log_with_order(
                &32942478.into(),
                &(Integer::from(11).pow(127)),
                &11.into(),
                &int("2745206")
            )
            .unwrap(),
            127
        );
    }

    #[test]
    fn known_factors() {
        // A prime modulus, so its factorization is the modulus itself.
        let n = int("432751500361");
        let n_factors = HashMap::from([(n.clone(), 1)]);
        assert_eq!(
            discrete_log_with_factors(&n, &(Integer::from(7).pow(324)), &7.into(), &n_factors)
                .unwrap(),
            324
        );
    }

    #[test]
    fn bases_sharing_a_factor_with_the_modulus() {
        // The powers of `b` are not a subgroup of the units, but the logarithm can still exist.
        assert_eq!(discrete_log(&9.into(), &0.into(), &3.into()).unwrap(), 2);
        assert_eq!(discrete_log(&12.into(), &4.into(), &2.into()).unwrap(), 2);
        assert_eq!(
            discrete_log(&10.into(), &3.into(), &2.into()),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn modulus_validation() {
        // Modulo 1 every number is 0.
        assert_eq!(discrete_log(&1.into(), &0.into(), &0.into()).unwrap(), 0);
        assert_eq!(
            discrete_log_with_order(&1.into(), &0.into(), &0.into(), &10.into()).unwrap(),
            0
        );

        for n in [0, -1, -4] {
            assert_eq!(
                discrete_log(&n.into(), &1.into(), &3.into()),
                Err(Error::InvalidModulus)
            );
            assert_eq!(
                discrete_log_with_order(&n.into(), &1.into(), &3.into(), &10.into()),
                Err(Error::InvalidModulus)
            );
            assert_eq!(
                discrete_log_with_prime_order(&n.into(), &1.into(), &3.into(), &11.into()),
                Err(Error::InvalidModulus)
            );
            assert_eq!(
                discrete_log_with_factors(&n.into(), &1.into(), &3.into(), &HashMap::new()),
                Err(Error::InvalidModulus)
            );
        }
    }
}
