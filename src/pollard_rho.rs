use rug::{Integer, rand::RandState};

use crate::{
    Error, check_modulus, discrete_log_trial_mul,
    modular::{BigRing, ModRing, WordRing},
    n_order,
};

/// Number of random starting points tried before giving up.
const RETRIES: usize = 10;

/// Number of multipliers of the walk.
///
/// An r-adding walk with that many multipliers is, in Teske's measurements, within a few percent
/// of a random mapping: fewer makes the walk collide later, more only costs the powers that build
/// them.
const MULTIPLIERS: usize = 20;

/// Multiplier of the partition: the odd number closest to `2**64 / phi`, whose multiples spread
/// over the whole word, the high bits of the product being the mixed ones.
const PARTITION_MULTIPLIER: u64 = 0x9E37_79B9_7F4A_7C15;

/// Steps of a walk per expected step of it, before it is abandoned for another one.
///
/// A walk over a group of `m` points closes a cycle after about `1.25 * sqrt(m)` steps, and Brent's
/// detection needs a few times that in the worst case: a walk still open long after that is
/// unlucky, and a fresh one is more likely to collide than its continuation.
const STEPS_PER_EXPECTED_STEP: u64 = 8;

/// One multiplier of the walk: `b**b_exponent * a**a_exponent`, and the exponents it adds.
struct Multiplier<G: ModRing, E: ModRing> {
    /// The group element the point is multiplied by.
    element: G::Elem,
    /// What the exponent of `b` grows by, modulo the order.
    b_exponent: E::Elem,
    /// What the exponent of `a` grows by, modulo the order.
    a_exponent: E::Elem,
}

/// Pollard's Rho algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// It is a randomized algorithm with the same expected running time as `discrete_log_shanks_steps`, but requires a negligible amount of memory.
///
/// The walk is an r-adding walk: the group is partitioned by the value of a point, and a point is
/// multiplied by the multiplier of its part, one of twenty random `b**u * a**v`. Its cycles are
/// found with Brent's algorithm, which walks a single point, one multiplication per step.
///
/// A modulus that is not positive is refused with [`Error::InvalidModulus`]. Modulo 1 every residue
/// is 0, so the logarithm is 0.
///
/// [`Error::LogDoesNotExist`] from here means "no logarithm found within the budget of the search",
/// not a proof that none exists: ten walks of a bounded number of steps each are tried, and a
/// logarithm that none of them met comes back as that error. Every candidate is verified before it
/// is returned, so a result is never wrong; only a failure can be.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation. It
/// must be the order of `b` modulo `n` for the result to be the smallest exponent: with a proper
/// multiple of it the walk returns whichever solution its relation yields, which is a valid
/// exponent and usually not the smallest one (`n = 7`, `b = 2`, `a = 4`, `order = 9` gives 8, not
/// 2). Sympy behaves the same way.
pub fn discrete_log_pollard_rho(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    solve(
        n,
        a,
        b,
        order,
        SequentialWalk {
            rand_state: &mut RandState::new(),
        },
    )
}

/// Pollard's Rho algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// Same as [`discrete_log_pollard_rho`], with the random generator seeded with `seed`: the same
/// seed always tries the same walks, and a retry with another seed makes other choices.
///
/// The preconditions and the meaning of every error are the ones of [`discrete_log_pollard_rho`].
pub fn discrete_log_pollard_rho_with_seed(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    seed: u64,
) -> Result<Integer, Error> {
    let mut rand_state = RandState::new();
    rand_state.seed(&Integer::from(seed));
    solve(
        n,
        a,
        b,
        order,
        SequentialWalk {
            rand_state: &mut rand_state,
        },
    )
}

/// What [`solve`] runs once it has chosen the rings the walk computes in.
///
/// The rings are chosen at run time while the walk is generic over them, so the choice cannot be
/// returned: it is the walk that is handed to the rings. The `Send` and `Sync` bounds are what the
/// parallel search needs of them, and both rings have them; the sequential walk ignores them.
trait RingTask {
    /// Solves the problem over the group `group`, the exponents being counted in `exponents`.
    ///
    /// `a` and `b` are already reduced modulo the modulus of `group`, and `order` is the order of
    /// `b`, at least four.
    fn run<G, E>(
        self,
        group: &G,
        exponents: &E,
        a: &Integer,
        b: &Integer,
        order: &Integer,
    ) -> Result<Integer, Error>
    where
        G: ModRing + Sync,
        G::Elem: Send + Sync,
        E: ModRing + Sync,
        E::Elem: Send + Sync;
}

/// Checks the inputs, finds the order of `b` when it is not given, and runs `task` over the rings
/// best suited to the modulus and to that order.
fn solve<T: RingTask>(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    task: T,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    // Modulo 1 every residue is 0, so every exponent is a logarithm and 0 is the smallest one. The
    // walk below would otherwise return whichever exponent its first relation happens to give.
    if *n == 1 {
        return Ok(Integer::new());
    }
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

    // The exponents are added and subtracted modulo the order at every step: a word-size ring for
    // them whenever the order fits in one, and then for the group too when the modulus does.
    match WordRing::from_modulus(&order) {
        Some(exponents) => match WordRing::from_modulus(n) {
            Some(group) => task.run(&group, &exponents, &a, &b, &order),
            None => {
                let group = BigRing::new(n).expect("a positive modulus");
                task.run(&group, &exponents, &a, &b, &order)
            }
        },
        None => {
            // An order above `2**64` is out of reach of a square root time algorithm anyway.
            let group = BigRing::new(n).expect("a positive modulus");
            let exponents = BigRing::new(&order).expect("an order of at least four");
            task.run(&group, &exponents, &a, &b, &order)
        }
    }
}

/// The task of [`discrete_log_pollard_rho`]: a single walk, in the calling thread.
struct SequentialWalk<'a, 'b> {
    /// Where the multipliers and the starting points are drawn from.
    rand_state: &'a mut RandState<'b>,
}

impl RingTask for SequentialWalk<'_, '_> {
    fn run<G, E>(
        self,
        group: &G,
        exponents: &E,
        a: &Integer,
        b: &Integer,
        order: &Integer,
    ) -> Result<Integer, Error>
    where
        G: ModRing + Sync,
        G::Elem: Send + Sync,
        E: ModRing + Sync,
        E::Elem: Send + Sync,
    {
        walk(group, exponents, a, b, order, self.rand_state)
    }
}

/// [`discrete_log_pollard_rho`] with the group in `group` and the exponents in `exponents`, the
/// order being at least four.
fn walk<G: ModRing, E: ModRing>(
    group: &G,
    exponents: &E,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error> {
    let a = group.from_integer(a);
    let b = group.from_integer(b);
    let budget = step_budget(order);

    for _ in 0..RETRIES {
        let multipliers = random_multipliers(group, exponents, &a, &b, order, rand_state);
        let (mut point, mut b_exponent, mut a_exponent) =
            random_start(group, exponents, &a, &b, order, rand_state);

        // Brent's cycle detection: a single point walks, compared at each step with the one saved
        // at the last power of two. Floyd's needs a second point walking twice as fast, three
        // multiplications per step instead of one.
        let mut saved = point.clone();
        let mut saved_b = b_exponent.clone();
        let mut saved_a = a_exponent.clone();
        let mut power = 1u64;
        let mut length = 0u64;

        for _ in 0..budget {
            advance(
                group,
                exponents,
                &multipliers,
                &mut point,
                &mut b_exponent,
                &mut a_exponent,
            );
            length += 1;

            if point == saved {
                if let Some(candidate) = relation(
                    group,
                    exponents,
                    &a,
                    &b,
                    &b_exponent,
                    &a_exponent,
                    &saved_b,
                    &saved_a,
                ) {
                    return Ok(candidate);
                }
                // The cycle holds no usable relation: only another walk can find one.
                break;
            }

            if length == power {
                // `clone_from` writes into the buffer the saved point already holds, where a
                // fresh clone would allocate one and free that one.
                saved.clone_from(&point);
                saved_b.clone_from(&b_exponent);
                saved_a.clone_from(&a_exponent);
                power *= 2;
                length = 0;
            }
        }
    }

    Err(Error::LogDoesNotExist)
}

/// How many steps a walk is given before it is abandoned, from the order of the group.
fn step_budget(order: &Integer) -> u64 {
    order
        .clone()
        .sqrt()
        .to_u64()
        .unwrap_or(u64::MAX)
        .saturating_mul(STEPS_PER_EXPECTED_STEP)
        .saturating_add(64)
}

/// The [`MULTIPLIERS`] multipliers of one walk, each a random `b**u * a**v`.
#[inline(never)]
fn random_multipliers<G: ModRing, E: ModRing>(
    group: &G,
    exponents: &E,
    a: &G::Elem,
    b: &G::Elem,
    order: &Integer,
    rand_state: &mut RandState<'_>,
) -> Vec<Multiplier<G, E>> {
    (0..MULTIPLIERS)
        .map(|_| {
            let (u, v) = random_exponents(order, rand_state);
            Multiplier {
                element: group.mul(&group.pow(b, &u), &group.pow(a, &v)),
                b_exponent: exponents.from_integer(&u),
                a_exponent: exponents.from_integer(&v),
            }
        })
        .collect()
}

/// A random point `b**u * a**v` to start a walk from, and the exponents it stands for.
#[inline(never)]
fn random_start<G: ModRing, E: ModRing>(
    group: &G,
    exponents: &E,
    a: &G::Elem,
    b: &G::Elem,
    order: &Integer,
    rand_state: &mut RandState<'_>,
) -> (G::Elem, E::Elem, E::Elem) {
    let (u, v) = random_exponents(order, rand_state);
    (
        group.mul(&group.pow(b, &u), &group.pow(a, &v)),
        exponents.from_integer(&u),
        exponents.from_integer(&v),
    )
}

/// A pair of exponents drawn below `order`.
fn random_exponents(order: &Integer, rand_state: &mut RandState<'_>) -> (Integer, Integer) {
    let u = Integer::from(order.random_below_ref(rand_state));
    let v = Integer::from(order.random_below_ref(rand_state));
    (u, v)
}

/// One step of the walk: the point is multiplied by the multiplier of its part, and the exponents
/// grow by the ones that multiplier stands for.
#[inline]
fn advance<G: ModRing, E: ModRing>(
    group: &G,
    exponents: &E,
    multipliers: &[Multiplier<G, E>],
    point: &mut G::Elem,
    b_exponent: &mut E::Elem,
    a_exponent: &mut E::Elem,
) {
    let multiplier = &multipliers[part_of(group, point)];
    group.mul_assign(point, &multiplier.element);
    exponents.add_assign(b_exponent, &multiplier.b_exponent);
    exponents.add_assign(a_exponent, &multiplier.a_exponent);
}

/// The logarithm two walks that met agree on, `None` when their relation holds none.
///
/// The two points are `b**bp * a**ap` and `b**bs * a**as`, and `a = b**x`: they are equal when
/// `bp - bs = x * (as - ap)` modulo the order. A denominator that is not invertible, and a
/// candidate that `b` to the power of is not `a`, are both useless and both possible: the order
/// need not be prime, and `a` need not be a power of `b` at all.
#[allow(clippy::too_many_arguments)]
fn relation<G: ModRing, E: ModRing>(
    group: &G,
    exponents: &E,
    a: &G::Elem,
    b: &G::Elem,
    b_exponent: &E::Elem,
    a_exponent: &E::Elem,
    saved_b: &E::Elem,
    saved_a: &E::Elem,
) -> Option<Integer> {
    let numerator = exponents.sub(b_exponent, saved_b);
    let denominator = exponents.sub(saved_a, a_exponent);
    let inverse = exponents.invert(&denominator)?;
    let candidate = exponents.to_integer(&exponents.mul(&numerator, &inverse));
    (group.pow(b, &candidate) == *a).then_some(candidate)
}

/// The part of the group `point` belongs to, one of [`MULTIPLIERS`].
///
/// The partition follows the value of the point, not the way the ring represents it, and mixes it
/// first: the low bits of a residue can be constant over a subgroup (every element of one modulo a
/// power of two is odd), where the high bits of a multiplication are not.
#[inline]
fn part_of<G: ModRing>(group: &G, point: &G::Elem) -> usize {
    let mixed = group.residue_word(point).wrapping_mul(PARTITION_MULTIPLIER);
    // The high bits of the product scaled to the number of parts, no division needed.
    ((u128::from(mixed) * MULTIPLIERS as u128) >> 64) as usize
}

#[cfg(test)]
mod tests {
    use rug::{integer::IsPrime, ops::Pow};

    use super::*;

    /// A prime modulus above `2**64` whose multiplicative group has a subgroup of order `order`,
    /// and a generator of that subgroup.
    fn big_modulus_and_base(order: &Integer) -> (Integer, Integer) {
        let mut multiple = (Integer::from(1) << 64u32) / order + 1u32;
        let n = loop {
            let candidate = Integer::from(&multiple * order) + 1u32;
            if candidate.is_probably_prime(30) != IsPrime::No {
                break candidate;
            }
            multiple += 1;
        };
        assert!(n > u64::MAX);

        let b = Integer::from(3)
            .pow_mod(&Integer::from(&n - 1u32).div_exact(order), &n)
            .unwrap();
        assert_ne!(b, 1, "3 is a {order}th power modulo {n}");
        (n, b)
    }

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
        // 5 is a primitive root modulo 227 and `5**132 = 3**7 (mod 227)`: sympy's walk, and the
        // three-branch walk ported from it, gave up on that logarithm, the r-adding walk finds it.
        assert_eq!(
            discrete_log_pollard_rho(&227.into(), &(Integer::from(3).pow(7)), &5.into(), None)
                .unwrap(),
            132
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

    #[test]
    fn modulus_above_a_word() {
        // The moduli GMP is kept for: the walk runs on `BigRing` and the exponents on `WordRing`.
        // A prime modulus above `2**64` with a subgroup small enough for a square root time walk.
        let order = Integer::from(1009u32);
        let (n, b) = big_modulus_and_base(&order);
        let a = b.clone().pow_mod(&Integer::from(225), &n).unwrap();
        assert_eq!(
            discrete_log_pollard_rho_with_seed(&n, &a, &b, Some(&order), 7).unwrap(),
            225
        );
    }

    #[test]
    fn invalid_modulus() {
        for n in [0, -1, -4] {
            assert_eq!(
                discrete_log_pollard_rho(&n.into(), &1.into(), &3.into(), None),
                Err(Error::InvalidModulus)
            );
        }
    }

    #[test]
    fn modulus_of_one() {
        // Modulo 1 every residue is 0, so the logarithm is 0. An order large enough to leave room
        // for a walk used to make it return whichever exponent the first relation gave.
        for order in [None, Some(Integer::from(1)), Some(Integer::from(1000))] {
            assert_eq!(
                discrete_log_pollard_rho(&1.into(), &0.into(), &0.into(), order.as_ref()).unwrap(),
                0,
                "order {order:?}"
            );
            assert_eq!(
                discrete_log_pollard_rho_with_seed(
                    &1.into(),
                    &5.into(),
                    &7.into(),
                    order.as_ref(),
                    3
                )
                .unwrap(),
                0,
                "order {order:?}"
            );
        }
    }

    #[test]
    fn order_that_is_a_multiple_of_the_real_one() {
        // The order of 2 modulo 7 is 3, and the walk is given 9: it returns a valid exponent, which
        // is not the smallest one. Sympy answers the same way, and this is what the documented
        // precondition on `order` is about.
        let x =
            discrete_log_pollard_rho_with_seed(&7.into(), &4.into(), &2.into(), Some(&9.into()), 1)
                .unwrap();
        assert_eq!(Integer::from(2).pow_mod(&x, &7.into()).unwrap(), 4);
        assert!(x < 9);
    }

    #[test]
    fn walk_in_a_known_group() {
        // Every logarithm of a small group, the walk run directly with both rings of words.
        let n = 1019u64;
        let order = 509u64;
        let group = WordRing::new(n).unwrap();
        let exponents = WordRing::new(order).unwrap();
        let b = Integer::from(4);
        let order = Integer::from(order);
        for x in [0u32, 1, 2, 17, 508] {
            let a = b
                .clone()
                .pow_mod(&Integer::from(x), &Integer::from(n))
                .unwrap();
            let mut rand_state = RandState::new();
            rand_state.seed(&Integer::from(x));
            assert_eq!(
                walk(&group, &exponents, &a, &b, &order, &mut rand_state).unwrap(),
                x,
                "log of {a} in base 4 modulo 1019"
            );
        }
    }

    #[test]
    fn parts_of_a_group() {
        // Every part of the partition is met over the group, and no value of it is out of range.
        let group = WordRing::new(1019).unwrap();
        let mut parts = vec![0usize; MULTIPLIERS];
        for value in 0..1019u64 {
            parts[part_of(&group, &group.from_u64(value))] += 1;
        }
        assert!(parts.iter().all(|&met| met > 0), "unused part: {parts:?}");
        // A partition that follows the value of an element and not its Montgomery form: the same
        // element of two rings of the same modulus lands in the same part.
        let big = BigRing::new(&Integer::from(1019)).unwrap();
        for value in 0..1019u64 {
            assert_eq!(
                part_of(&group, &group.from_u64(value)),
                part_of(&big, &big.from_integer(&Integer::from(value)))
            );
        }
    }

    #[test]
    fn walk_exponents() {
        // A point of the walk always stays `b**u * a**v` of the exponents the walk adds up.
        let n = Integer::from(1019);
        let order = Integer::from(509);
        let group = WordRing::new(1019).unwrap();
        let exponents = WordRing::new(509).unwrap();
        let (a, b) = (Integer::from(625), Integer::from(4));
        let mut rand_state = RandState::new();
        rand_state.seed(&Integer::from(3));

        let a_elem = group.from_integer(&a);
        let b_elem = group.from_integer(&b);
        let multipliers = random_multipliers(
            &group,
            &exponents,
            &a_elem,
            &b_elem,
            &order,
            &mut rand_state,
        );

        let mut point = group.one();
        let mut b_exponent = exponents.zero();
        let mut a_exponent = exponents.zero();
        for _ in 0..2000 {
            advance(
                &group,
                &exponents,
                &multipliers,
                &mut point,
                &mut b_exponent,
                &mut a_exponent,
            );

            let expected = b
                .clone()
                .pow_mod(&exponents.to_integer(&b_exponent), &n)
                .unwrap()
                * a.clone()
                    .pow_mod(&exponents.to_integer(&a_exponent), &n)
                    .unwrap()
                % &n;
            assert_eq!(group.to_integer(&point), expected);
        }
    }
}
