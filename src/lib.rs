#![doc = include_str!("../README.md")]
//!
//! # Optional features
//!
//! - `parallel` (off by default): adds `discrete_log_pollard_rho_parallel`, the parallel collision
//!   search of van Oorschot and Wiener, whose expected running time is that of Pollard's Rho
//!   divided by the number of threads. It is opt-in because it is the only entry point that spawns
//!   threads, [`discrete_log`] and the others staying in the thread that calls them. It needs no
//!   dependency, `std::thread` only.
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use std::collections::HashMap;

use rug::{Integer, integer::IsPrime, ops::Pow};

use factor::fast_factor;
use n_order::element_order_with_factors;

mod crt;
mod factor;
mod index_calculus;
mod modular;
mod n_order;
mod pohlig_hellman;
mod pollard_rho;
mod shanks_steps;
mod trial_mul;

pub use index_calculus::{discrete_log_index_calculus, discrete_log_index_calculus_with_seed};
pub use n_order::{n_order, n_order_with_factors};
pub use pohlig_hellman::{discrete_log_pohlig_hellman, discrete_log_pohlig_hellman_with_factors};
#[cfg(feature = "parallel")]
pub use pollard_rho::{MAX_THREADS, discrete_log_pollard_rho_parallel};
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
///
/// The type is `#[non_exhaustive]`: a later version can tell a new failure apart from the ones
/// below without that being a breaking change, so a match on it needs a catch-all arm.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// Log does not exist
    ///
    /// From the randomized algorithms it means "no logarithm found within the budget of the
    /// search", which is not a proof that none exists: see [`discrete_log_pollard_rho`] and
    /// [`discrete_log_index_calculus`].
    #[error("Log does not exist")]
    LogDoesNotExist,
    /// A and n are not relatively prime
    #[error("A and n are not relatively prime")]
    NotRelativelyPrime,
    /// The modulus is not positive
    #[error("n should be positive")]
    InvalidModulus,
    /// The order given, or the factorization given of it, cannot be one
    ///
    /// An order that is not positive, or a factorization holding a prime below 2, an exponent below
    /// 1, or prime powers that do not multiply back to the number it is supposed to factor. The
    /// entry points taking a factorization ([`discrete_log_with_factors`],
    /// [`discrete_log_pohlig_hellman_with_factors`], [`n_order_with_factors`]) check it before they
    /// compute anything: what follows would otherwise divide by zero, raise to a negative exponent
    /// or loop without end.
    #[error("the order or its factorization is not usable")]
    InvalidOrder,
    /// The instance is beyond what the chosen algorithm can do
    ///
    /// Unlike [`Error::LogDoesNotExist`] it says nothing about the logarithm: it may well exist,
    /// and another algorithm may well find it. [`discrete_log_shanks_steps`] returns it for an
    /// order whose table of baby steps is above its memory cap, and
    /// [`discrete_log_index_calculus`] for a modulus whose relation matrix is above its own.
    #[error("the instance is out of reach of this algorithm")]
    OutOfReach,
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
    // The factorization is this crate's own: nothing to check about it.
    solve_with_factors(n, a, b, &fast_factor(n))
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// If the prime factorization of `n` is known, it can be passed as `n_factors` to speed up the computation.
///
/// `n_factors` must be the prime factorization of `n`, its primes mapped to their exponents.
/// Anything else is refused with [`Error::InvalidOrder`], the prime powers being multiplied back
/// together first: the order of the base is derived from that map, and a map belonging to another
/// number gives no order at all.
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
    check_factorization(n, n_factors)?;
    solve_with_factors(n, a, b, n_factors)
}

/// [`discrete_log_with_factors`] with `n_factors` already known to be the factorization of `n`, and
/// `n` already known to be above 1.
fn solve_with_factors(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    n_factors: &HashMap<Integer, usize>,
) -> Result<Integer, Error> {
    let b = b.clone().modulo(n);
    let (order, order_factors) = element_order_with_factors(n, &b, n_factors);
    // A single prime with exponent 1: the order is prime, no primality test needed.
    let prime_order = order_factors.len() == 1 && order_factors.values().all(|&e| e == 1);
    solve(n, a, &b, &order, Some(&order_factors), Some(prime_order))
}

/// Compute the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation.
///
/// `order` must be the order of `b` modulo `n` for the result to be the smallest exponent. Any
/// multiple of it still leads to an exponent `x` with `b**x = a (mod n)`, since the algorithms
/// search `0..order`, but which one depends on the algorithm the order sends the problem to:
/// baby-step giant-step returns the smallest solution whatever the multiple is, Pollard's rho
/// returns whichever solution its walk yields, and Pohlig-Hellman can return
/// [`Error::LogDoesNotExist`] for a logarithm that exists (`n = 7`, `b = 2`, `a = 4`, `order = 9`).
/// Sympy behaves the same way.
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
///
/// The precondition on `order` is the one of [`discrete_log_with_order`]: only the real order of
/// `b` guarantees the smallest exponent.
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

/// Orders below this one are searched exhaustively, one power of the base at a time.
///
/// Sympy's threshold, kept: every other algorithm builds a table, a walk or a factor base before it
/// looks at a single exponent, which is already most of what walking through a few hundred
/// exponents costs. Nothing here measures trial multiplication, so there is nothing to retune it
/// from either.
pub(crate) const TRIAL_MUL_ORDER: u64 = 1000;

/// Prime orders below this one go to baby-step giant-step, the larger ones to Pollard's rho.
///
/// The two are both `O(sqrt(order))`, so only their constants and their memory decide. Measured
/// over twelve safe-prime instances per size (three logarithms each), the two meet at an order of
/// `2**31`: in estimated cycles baby-step giant-step is 1.3 times cheaper at `2**30` and 1.2 times
/// dearer at `2**33`. Its instruction count alone stays the lower one up to about `2**34`, but the
/// table it walks is what its instructions wait for: its last-level cache misses grow with the
/// table (3.1M at an order of `2**29`, 18.5M at `2**34`) where Pollard's rho keeps about 20k
/// whatever the order is, and wall-clock time, which pays for those misses in full, already
/// favours rho from `2**29` on. `2**31` is where the three measures agree the most, and it is also
/// where the table is still a megabyte or so, against no memory at all for rho.
///
/// Below [`shanks_steps::MAX_ORDER`], the size the table stops fitting in memory at, so that the
/// larger orders are never handed to an algorithm that refuses them.
pub(crate) const SHANKS_STEPS_ORDER: u64 = 1 << 31;

const _: () = assert!(SHANKS_STEPS_ORDER <= shanks_steps::MAX_ORDER);

/// Slack of the index calculus comparison, in natural logarithms of the group order.
///
/// Index calculus is subexponential in the modulus where the two square root algorithms are
/// exponential in the order, so the costs to compare are `exp(2*sqrt(log(n)*log(log(n))))` against
/// `sqrt(order)`, whose logarithms are `2*sqrt(log(n)*log(log(n)))` and `log(order)/2`. Neither
/// carries its implementation constant, and this is what stands in for the ratio of the two:
/// doubling both sides, index calculus is chosen when `4*sqrt(log(n)*log(log(n)))` is below
/// `log(order)` plus this.
///
/// Sympy uses `-10`, which on a safe prime puts the boundary above 130 bits. Measured over safe
/// primes of 34 to 56 bits, three logarithms each, index calculus overtakes Pollard's rho at about
/// 35 bits, where this expression is worth 11.6, and is twice as cheap from 44 bits on, where it is
/// worth 11.0. The measured crossover is not what this is set to: index calculus can give up where
/// rho cannot (a lost relation search is a [`Error::LogDoesNotExist`] for a logarithm that exists),
/// so it is only worth choosing where it wins by a clear factor. 11.0 is that boundary.
///
/// One unit of it is worth seven to eight bits of order, not five: along the safe prime family,
/// where the order is half the modulus, the boundary falls at an order of about 60 bits with a
/// slack of 9, 53 bits with 10 and about 45 bits with 11.
pub(crate) const INDEX_CALCULUS_SLACK: f64 = 11.0;

/// Whether index calculus is expected to beat the square root algorithms on a prime `order`
/// modulo `n`, by the comparison [`INDEX_CALCULUS_SLACK`] describes.
///
/// A modulus too large for a `f64` makes the left hand side infinite and so the answer `false`,
/// which is the only useful one: the factor base of such a modulus is far beyond what can be
/// sieved, and [`discrete_log_index_calculus`] refuses it.
pub(crate) fn index_calculus_pays(n: &Integer, order: &Integer) -> bool {
    let log_n = n.to_f64().ln();
    let log_order = order.to_f64().ln();
    4.0 * (log_n * log_n.ln()).sqrt() < log_order + INDEX_CALCULUS_SLACK
}

/// The algorithms [`solve`] picks from for a base of prime order.
pub(crate) enum Algorithm {
    /// Exhaustive search over the powers of the base.
    TrialMul,
    /// Index calculus, subexponential in the modulus.
    IndexCalculus,
    /// Baby-step giant-step, `O(sqrt(order))` time and memory.
    ShanksSteps,
    /// Pollard's rho, `O(sqrt(order))` time and no memory.
    PollardRho,
}

/// The algorithm best suited to a base of prime `order` modulo `n`.
///
/// The one place the choice is made: anything that stands in for one of these algorithms, like the
/// shared table of baby steps of Pohlig-Hellman, asks here instead of comparing sizes of its own.
///
/// The bounds on the order are tested before [`index_calculus_pays`], and not after it as sympy
/// does. That comparison is a difference of two growth rates with no implementation constant in it,
/// and `4*sqrt(log(n)*log(log(n)))` peaks in slope around `log(n) = 20`: with the slack needed for
/// the orders it is meant to decide, it also comes out true for the small moduli, where measurement
/// puts index calculus 3.6 to 6.6 times behind the baby-step giant-step it would displace (8.94 us
/// against 1.36 us at `n = 2027`, 34.14 us against 5.37 us at `n = 100043`). An order below
/// [`SHANKS_STEPS_ORDER`] is a table of at most a megabyte and a search of at most `2**16` steps,
/// which nothing subexponential in the modulus beats, so it never reaches the comparison.
pub(crate) fn algorithm_for(n: &Integer, order: &Integer) -> Algorithm {
    if *order < TRIAL_MUL_ORDER {
        Algorithm::TrialMul
    } else if *order < SHANKS_STEPS_ORDER {
        Algorithm::ShanksSteps
    } else if index_calculus_pays(n, order) {
        Algorithm::IndexCalculus
    } else {
        Algorithm::PollardRho
    }
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
    if *order < TRIAL_MUL_ORDER {
        return discrete_log_trial_mul(n, a, b, Some(order));
    }

    // GMP runs a Baillie-PSW test and this many Miller-Rabin rounds on top of it, each of which a
    // composite passes with probability at most 1/4: the factorization uses the same number, and
    // nothing here is fast enough for the handful of modular powerings they cost to show.
    let prime_order = prime_order
        .unwrap_or_else(|| order.is_probably_prime(factor::PRIMALITY_REPS) != IsPrime::No);
    if !prime_order {
        return match order_factors {
            // The factorization comes from `element_order_with_factors`: nothing to check about it.
            Some(order_factors) => {
                pohlig_hellman::solve_with_factors(n, a, b, order, order_factors)
            }
            None => discrete_log_pohlig_hellman(n, a, b, Some(order)),
        };
    }

    match algorithm_for(n, order) {
        // The exhaustive search above, which a prime order below the threshold has already taken.
        Algorithm::TrialMul => discrete_log_trial_mul(n, a, b, Some(order)),
        Algorithm::IndexCalculus => discrete_log_index_calculus(n, a, b, Some(order)),
        Algorithm::ShanksSteps => discrete_log_shanks_steps(n, a, b, Some(order)),
        Algorithm::PollardRho => discrete_log_pollard_rho(n, a, b, Some(order)),
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

/// Rejects a map that cannot be the prime factorization of `n`.
///
/// The primes must be at least 2, their exponents at least 1, and the prime powers must multiply
/// back to `n`, which must itself be positive. That is one powering and one multiplication per
/// factor, nothing next to the algorithm it guards, and it is what keeps the factorization of
/// another number out: the code that trusts such a map divides by zero, raises to a negative
/// exponent, or walks a loop that never ends because the base it was given is not a unit of the
/// prime power it is asked about.
///
/// Primality of the primes is the one part of the contract left to the caller: a Miller-Rabin
/// round per factor costs more than everything checked here, and a composite key still multiplies
/// out to `n` as often as not.
pub(crate) fn check_factorization(
    n: &Integer,
    factors: &HashMap<Integer, usize>,
) -> Result<(), Error> {
    if *n < 1 {
        return Err(Error::InvalidOrder);
    }

    let mut product = Integer::from(1);
    for (p, e) in factors {
        // A prime of at least 2 raised to the number of bits of `n` is already above `n`, so the
        // product could never come back to it: refusing the exponent here is also what keeps the
        // cast below from wrapping one above `u32::MAX` down to a plausible power.
        if *p < 2 || *e < 1 || *e as u64 >= u64::from(n.significant_bits()) {
            return Err(Error::InvalidOrder);
        }
        product *= p.clone().pow(*e as u32);
        // The factors are multiplied in whatever order the map holds them, and the product only
        // grows: one that is already too large is one that will stay too large.
        if product > *n {
            return Err(Error::InvalidOrder);
        }
    }

    if product == *n {
        Ok(())
    } else {
        Err(Error::InvalidOrder)
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
        let n = int(
            "83408372012221120677052349409462320990177094246143674474872152829440524098582262384066400107950985845255268335597502228206679771838750219696329523257176739436871327238322817403970284015587320158034304282786944710043150568360761457471641695390427267786485448748458445872307883254297662715749746270343116946519",
        );
        let a = int(
            "109770827223661560471527567179288748906402603483328748683689436879660543465776899146036833470531024202351087008847594392666852763100570391337823820240726499421306887565697452868723849092658743267256316770223643723095601213088336064635680075206929620159782416078143076506249031972043819429093074684182845530529249907297736582589125917235222921623698038868900282049587768700860009877737045693722732170123306528145661683416808514556360429554775212088169626620488741903267154641722293484797745665402402381445609873333905772582972140944493849645600529147490903067975300304532955461710562911203871840101407995813072692212",
        );
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
    fn factorization_validation() {
        // What is and is not the prime factorization of a number.
        for (n, factors, ok) in [
            (1, vec![], true),
            (9, vec![(3, 2)], true),
            (12, vec![(2, 2), (3, 1)], true),
            // A composite key that still multiplies back: the primality of the keys is the part of
            // the contract this does not check.
            (4, vec![(4, 1)], true),
            // Nothing is a factorization of a number that is not positive.
            (0, vec![], false),
            (-4, vec![(2, 2)], false),
            // Primes below 2, and negative ones.
            (9, vec![(0, 1)], false),
            (9, vec![(1, 1)], false),
            (9, vec![(-3, 2)], false),
            // Exponents below 1, and one that a `u32` would wrap.
            (9, vec![(3, 0)], false),
            (9, vec![(3, usize::MAX)], false),
            // The prime powers of another number.
            (9, vec![(3, 1)], false),
            (9, vec![(3, 3)], false),
            (9, vec![(3, 2), (5, 1)], false),
            (1, vec![(2, 1)], false),
        ] {
            let map: HashMap<Integer, usize> = factors
                .iter()
                .map(|&(p, e)| (Integer::from(p), e))
                .collect();
            assert_eq!(
                check_factorization(&Integer::from(n), &map).is_ok(),
                ok,
                "n = {n}, factors = {factors:?}"
            );
        }

        // And how a wrong one reaches the caller.
        assert_eq!(
            discrete_log_with_factors(
                &9.into(),
                &1.into(),
                &2.into(),
                &HashMap::from([(Integer::from(3), 0)])
            ),
            Err(Error::InvalidOrder)
        );
        assert_eq!(
            discrete_log_with_factors(
                &12.into(),
                &4.into(),
                &2.into(),
                &HashMap::from([(Integer::from(2), 0)])
            ),
            Err(Error::InvalidOrder)
        );
        // The factorization that does belong to the modulus still solves.
        assert_eq!(
            discrete_log_with_factors(
                &12.into(),
                &4.into(),
                &2.into(),
                &HashMap::from([(Integer::from(2), 2), (Integer::from(3), 1)])
            )
            .unwrap(),
            2
        );
    }

    #[test]
    fn algorithm_selection() {
        // A prime order below `SHANKS_STEPS_ORDER` goes to baby-step giant-step, although the index
        // calculus comparison comes out true for these moduli: it is a table of at most a megabyte
        // and a search of at most `2**16` steps, which measurement puts 3.6 to 6.6 times ahead of
        // the index calculus that used to displace it.
        for (n, order) in [(2027u32, 1013u32), (100_043, 50_021), (1_158_197, 579_098)] {
            let (n, order) = (Integer::from(n), Integer::from(order));
            assert!(index_calculus_pays(&n, &order), "n = {n}");
            assert!(
                matches!(algorithm_for(&n, &order), Algorithm::ShanksSteps),
                "n = {n}"
            );
        }

        // Below the exhaustive search threshold nothing else is asked.
        assert!(matches!(
            algorithm_for(&2027.into(), &999.into()),
            Algorithm::TrialMul
        ));
        // Above the table, the comparison decides: a safe prime of 72 bits is for index calculus, a
        // 41-bit order under a 41-bit modulus for Pollard's rho.
        assert!(matches!(
            algorithm_for(
                &int("4169048128772662588679"),
                &int("2084524064386331294339")
            ),
            Algorithm::IndexCalculus
        ));
        assert!(matches!(
            algorithm_for(&int("2199023255867"), &int("1099511627933")),
            Algorithm::PollardRho
        ));
    }

    #[test]
    fn modulus_out_of_reach() {
        // A safe prime of 128 bits: the order is prime and far above the table of baby steps, so the
        // dispatcher sends it to index calculus, whose factor base would need 2.4 GiB of relation
        // matrix. It used to spin over candidates for that matrix and never return.
        let n = int("279305127720133152706028210119365912959");
        let a = int("269390257626015088577421517053028557938");
        assert_eq!(discrete_log(&n, &a, &4.into()), Err(Error::OutOfReach));
        assert_eq!(
            discrete_log_with_prime_order(&n, &a, &4.into(), &(Integer::from(&n - 1u32) / 2u32)),
            Err(Error::OutOfReach)
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
