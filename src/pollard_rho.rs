#[cfg(feature = "parallel")]
use std::{
    collections::HashMap,
    sync::{
        Mutex, MutexGuard,
        atomic::{AtomicBool, Ordering},
    },
};

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

/// Multiplier of the distinguished point test: another odd number with its bits well spread, and
/// not the one the partition uses, so that being distinguished tells nothing about the part a
/// point falls in.
#[cfg(feature = "parallel")]
const DISTINGUISHED_MULTIPLIER: u64 = 0xC2B2_AE3D_27D4_EB4F;

/// Most leading zeros a distinguished point is asked for, whatever the order.
///
/// One point in `2**24` is distinguished, so a walk covers sixteen million steps between two of
/// them: past that the table stops growing usefully and the tail each thread walks after the
/// collision starts to be felt.
#[cfg(feature = "parallel")]
const MAX_DISTINGUISHED_BITS: u32 = 24;

/// Most distinguished points the shared table holds.
///
/// A point costs its residue and two exponents: three words for a modulus that fits in one, a
/// dozen megabytes of table once the load factor is counted, and three GMP integers otherwise, a
/// few tens of megabytes. With [`distinguished_bits`] the expected number of points stored is the
/// fourth root of the order, below the cap for every order a square root time search can reach, so
/// only a search that finds nothing and walks its whole budget fills the table.
#[cfg(feature = "parallel")]
const MAX_DISTINGUISHED_POINTS: usize = 1 << 18;

/// Steps between two reads of the stop signal by a walk that meets no distinguished point.
#[cfg(feature = "parallel")]
const STOP_CHECK_STEPS: u64 = 1024;

/// Most walks a parallel search runs, whatever it is asked for.
///
/// A count taken as it comes is a count the caller can make the process die of: a hundred thousand
/// threads is a spawn the operating system refuses (`WouldBlock`, a panic out of
/// `std::thread::scope`), and `usize::MAX` of them overflows the vector of seeds before a single
/// walk starts. No machine this runs on has this many cores, and more walks than there are cores
/// only makes each of them slower, so clamping here costs nothing that was ever worth having.
#[cfg(feature = "parallel")]
pub const MAX_THREADS: usize = 256;

/// The exponents of `b` and of `a` a point was reached with.
#[cfg(feature = "parallel")]
type Exponents<E> = (<E as ModRing>::Elem, <E as ModRing>::Elem);

/// Pollard's Rho algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`), searched by `threads` threads at once.
///
/// This is the parallel collision search of van Oorschot and Wiener, and the only entry point of
/// the crate that spawns threads: every other one stays in the thread that calls it. It is opt-in
/// for that reason, behind the `parallel` feature.
///
/// The threads walk the same r-adding walk as [`discrete_log_pollard_rho`], each from a starting
/// point of its own. A point is *distinguished* when a cheap hash of its residue comes out with
/// enough leading zeros, and the distinguished points go into a table shared by the threads,
/// together with the exponents they were reached with. Two walks that meet are the same walk from
/// there on, so they reach the same distinguished point, and the two exponent pairs it was stored
/// with give the logarithm. Waiting for a walk to close its own cycle instead, as the sequential
/// version does, gives each thread the work of a whole search: this divides the expected number of
/// steps by the number of threads, which no sequential variant of the walk can do.
///
/// One point in `order**(1/4)` is distinguished, so the table holds about `order**(1/4)` points
/// against the `order**(1/2)` steps of the search, and each thread walks about `order**(1/4)` steps
/// past the collision before it notices. It is capped anyway, at a few hundred thousand points, a
/// few tens of megabytes at worst.
///
/// `threads` of 0 or 1 runs [`discrete_log_pollard_rho`] in the calling thread instead. More
/// threads than the machine has cores only slows the search down, and an order small enough for
/// the walk to end in microseconds is not worth spawning a thread for. A count above
/// [`MAX_THREADS`] is clamped to it rather than refused, and a spawn the operating system will not
/// grant leaves the search to the walks already started: asking for a hundred thousand threads is
/// answered by a search over [`MAX_THREADS`] of them, not by a panic.
///
/// One of the walks runs in the calling thread, which would otherwise only wait for the others, so
/// `threads` walks cost `threads - 1` spawned threads.
///
/// A modulus that is not positive is refused with [`Error::InvalidModulus`]. Modulo 1 every residue
/// is 0, so the logarithm is 0. [`Error::LogDoesNotExist`] means "no logarithm found within the
/// budget of the search", as it does for [`discrete_log_pollard_rho`]: every thread walks a bounded
/// number of steps, so a problem with no logarithm ends in a failure instead of a search that never
/// stops.
///
/// Threads make the walks unreproducible, so which relation solves the problem changes from one run
/// to the next. The logarithm is verified before it is returned, and reduced into `[0, order)`.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation. The
/// precondition on it is the one of [`discrete_log_pollard_rho`]: only the real order of `b`
/// guarantees the smallest exponent.
///
/// # Examples
///
/// ```
/// use discrete_logarithm::discrete_log_pollard_rho_parallel;
/// use rug::Integer;
///
/// // A prime modulus whose multiplicative group has a subgroup of prime order `2**30`.
/// let n = Integer::from(2147483783u32);
/// let order = Integer::from(1073741891u32);
/// let b = Integer::from(9);
/// let a = Integer::from(240127149u32);
/// let x = discrete_log_pollard_rho_parallel(&n, &a, &b, Some(&order), 4).unwrap();
/// assert_eq!(x, 12345678);
/// ```
#[cfg(feature = "parallel")]
pub fn discrete_log_pollard_rho_parallel(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    threads: usize,
) -> Result<Integer, Error> {
    let mut rand_state = RandState::new();
    let threads = threads.min(MAX_THREADS);
    if threads <= 1 {
        return solve(
            n,
            a,
            b,
            order,
            SequentialWalk {
                rand_state: &mut rand_state,
            },
        );
    }
    solve(
        n,
        a,
        b,
        order,
        ParallelSearch {
            threads,
            rand_state: &mut rand_state,
        },
    )
}

/// The task of [`discrete_log_pollard_rho_parallel`]: several walks, one per spawned thread.
#[cfg(feature = "parallel")]
struct ParallelSearch<'a, 'b> {
    /// How many threads walk, at least two.
    threads: usize,
    /// Where the multipliers, the starting points and the seeds of the threads are drawn from.
    rand_state: &'a mut RandState<'b>,
}

#[cfg(feature = "parallel")]
impl RingTask for ParallelSearch<'_, '_> {
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
        parallel_walk(group, exponents, a, b, order, self.threads, self.rand_state)
    }
}

/// What the threads of a [`Search`] write to.
#[cfg(feature = "parallel")]
struct Shared<G: ModRing, E: ModRing> {
    /// Every distinguished point met so far, and the exponents the walk that stored it had there.
    table: Mutex<HashMap<G::Elem, Exponents<E>>>,
    /// The logarithm, once a thread has derived one from a collision and verified it.
    answer: Mutex<Option<Integer>>,
    /// Set once the search is over, so that the other threads leave their walk.
    stop: AtomicBool,
}

/// One parallel collision search: what every thread reads, and what they all write to.
#[cfg(feature = "parallel")]
struct Search<'a, G: ModRing, E: ModRing> {
    /// The group the walk multiplies in.
    group: &'a G,
    /// The ring the exponents are counted in, of modulus the order.
    exponents: &'a E,
    /// The element whose logarithm is wanted.
    a: G::Elem,
    /// The base of the logarithm.
    b: G::Elem,
    /// The order of `b`, at least four.
    order: &'a Integer,
    /// The multipliers of the walk, the same ones for every thread.
    multipliers: Vec<Multiplier<G, E>>,
    /// Steps one thread walks from one starting point before trying another.
    budget: u64,
    /// Leading zeros of the hash of a residue that make its point distinguished.
    distinguished_bits: u32,
    /// What the threads write to.
    shared: Shared<G, E>,
}

/// [`discrete_log_pollard_rho_parallel`] with the group in `group` and the exponents in
/// `exponents`, the order being at least four and `threads` being at least two.
#[cfg(feature = "parallel")]
fn parallel_walk<G, E>(
    group: &G,
    exponents: &E,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    threads: usize,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error>
where
    G: ModRing + Sync,
    G::Elem: Send + Sync,
    E: ModRing + Sync,
    E::Elem: Send + Sync,
{
    let a = group.from_integer(a);
    let b = group.from_integer(b);
    // One set of multipliers for the whole search, where the sequential walk draws fresh ones at
    // every retry: a collision between two walks is only a relation if both add the same
    // exponents, and new multipliers would make the points already in the table meaningless.
    let multipliers = random_multipliers(group, exponents, &a, &b, order, rand_state);
    // A generator cannot cross a thread boundary, so each thread builds its own from a seed drawn
    // here: the whole search still comes out of this one generator, and the walks stay independent.
    let seeds: Vec<u64> = (0..threads).map(|_| random_word(rand_state)).collect();

    let search = Search {
        group,
        exponents,
        a,
        b,
        order,
        multipliers,
        budget: step_budget(order),
        distinguished_bits: distinguished_bits(order),
        shared: Shared {
            table: Mutex::new(HashMap::new()),
            answer: Mutex::new(None),
            stop: AtomicBool::new(false),
        },
    };

    // A scope borrows the search instead of counting references to it, and joins every thread
    // before it returns: there is no thread left running when the table and the answer go away.
    std::thread::scope(|scope| {
        let search = &search;
        // The first seed is walked here, in the thread that would otherwise only wait for the
        // others, so `threads` walks need `threads - 1` spawns. A spawn the operating system will
        // not grant stops the loop instead of unwinding out of the scope: the walks already started
        // are then the whole search, which is a slower search and not a failed one.
        let (mine, spawned) = seeds.split_first().expect("at least two seeds");
        for &seed in spawned {
            if std::thread::Builder::new()
                .spawn_scoped(scope, move || search.trail(seed))
                .is_err()
            {
                break;
            }
        }
        search.trail(*mine);
    });

    match take(&search.shared.answer) {
        Some(candidate) => Ok(candidate),
        None => Err(Error::LogDoesNotExist),
    }
}

#[cfg(feature = "parallel")]
impl<G: ModRing, E: ModRing> Search<'_, G, E> {
    /// The walks of one thread, from starting points drawn off `seed`.
    ///
    /// It returns once it has derived the logarithm, once another thread has, or once every
    /// starting point it was given has run out of steps. Each of those is bounded, so a problem
    /// with no logarithm ends the thread instead of walking forever.
    fn trail(&self, seed: u64) {
        let mut rand_state = RandState::new();
        rand_state.seed(&Integer::from(seed));

        for _ in 0..RETRIES {
            if self.stopped() {
                return;
            }
            let (mut point, mut b_exponent, mut a_exponent) = random_start(
                self.group,
                self.exponents,
                &self.a,
                &self.b,
                self.order,
                &mut rand_state,
            );

            for step in 0..self.budget {
                advance(
                    self.group,
                    self.exponents,
                    &self.multipliers,
                    &mut point,
                    &mut b_exponent,
                    &mut a_exponent,
                );

                if !self.is_distinguished(&point) {
                    // The stop signal is read at the distinguished points, which is where the
                    // shared table is touched anyway; a walk can go a long way between two of
                    // them, and this bounds how long it keeps going after the search is over
                    // without an atomic read at every step.
                    if step % STOP_CHECK_STEPS == 0 && self.stopped() {
                        return;
                    }
                    continue;
                }

                let Some((saved_b, saved_a)) = self.store(&point, &b_exponent, &a_exponent) else {
                    continue;
                };
                // Derived and verified with the table unlocked: the exponentiation the check costs
                // is longer than the lookup, and every other thread would be waiting on it.
                if let Some(candidate) = relation(
                    self.group,
                    self.exponents,
                    &self.a,
                    &self.b,
                    &b_exponent,
                    &a_exponent,
                    &saved_b,
                    &saved_a,
                ) {
                    let mut answer = lock(&self.shared.answer);
                    // Two threads can finish at once, and either logarithm is verified: the first
                    // one stands, so that the answer does not depend on which lock was won.
                    if answer.is_none() {
                        *answer = Some(candidate);
                    }
                    drop(answer);
                    self.shared.stop.store(true, Ordering::Release);
                    return;
                }
                // A collision with no usable relation in it: the walk goes on to the next one.
                if self.stopped() {
                    return;
                }
            }
        }
    }

    /// Whether `point` is one of the distinguished points, the ones the threads share.
    ///
    /// The test is a multiplication and a shift, cheap enough to run at every step of the walk.
    /// The hash is of the residue and not of the representation, as [`part_of`] is, so that the
    /// same point is distinguished in whichever ring the group happens to be held.
    #[inline]
    fn is_distinguished(&self, point: &G::Elem) -> bool {
        let hash = self
            .group
            .residue_word(point)
            .wrapping_mul(DISTINGUISHED_MULTIPLIER);
        hash >> (u64::BITS - self.distinguished_bits) == 0
    }

    /// Stores a distinguished point and the exponents it was reached with, returning the exponents
    /// an earlier walk had at that same point when there is one.
    ///
    /// A full table stores nothing more: the search then only sees the collisions that involve a
    /// point already in it, which is all of them until the run is far longer than its expected
    /// length.
    fn store(
        &self,
        point: &G::Elem,
        b_exponent: &E::Elem,
        a_exponent: &E::Elem,
    ) -> Option<Exponents<E>> {
        let mut table = lock(&self.shared.table);
        match table.get(point) {
            Some(exponents) => Some(exponents.clone()),
            None => {
                if table.len() < MAX_DISTINGUISHED_POINTS {
                    table.insert(point.clone(), (b_exponent.clone(), a_exponent.clone()));
                }
                None
            }
        }
    }

    /// Whether the search is over, here or in another thread.
    #[inline]
    fn stopped(&self) -> bool {
        self.shared.stop.load(Ordering::Acquire)
    }
}

/// Leading zeros the hash of a residue must have for its point to be distinguished.
///
/// One point in `2**bits` is distinguished, which trades the memory of the table against the work
/// thrown away at the end of the search: a run of `s` steps stores about `s / 2**bits` points, and
/// each thread walks `2**bits` steps on average after the collision that solves the problem before
/// it reaches the distinguished point that shows it. A quarter of the bits of the order makes both
/// of them the fourth root of the order, against the square root of it the search itself costs.
#[cfg(feature = "parallel")]
fn distinguished_bits(order: &Integer) -> u32 {
    (order.significant_bits() / 4).clamp(1, MAX_DISTINGUISHED_BITS)
}

/// A word drawn from `rand_state`, to seed the generator of one thread with.
#[cfg(feature = "parallel")]
fn random_word(rand_state: &mut RandState<'_>) -> u64 {
    Integer::from(Integer::random_bits(u64::BITS, rand_state)).to_u64_wrapping()
}

/// Locks `mutex`, taking its value back when a panic elsewhere poisoned it.
///
/// Nothing in the search panics while holding a lock, so what is behind it is whole whatever the
/// poison says, and bringing down a thread of the search over it would only lose the walks the
/// others have done.
#[cfg(feature = "parallel")]
fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// Takes the value out of `mutex`, poisoned or not.
#[cfg(feature = "parallel")]
fn take<T: Default>(mutex: &Mutex<T>) -> T {
    std::mem::take(&mut *lock(mutex))
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

    /// Two instances of prime order just above `2**40`: `b` generates the subgroup of that order
    /// of the multiplicative group of a safe prime `2 * order + 1`.
    #[cfg(feature = "parallel")]
    const PRIME_ORDER_INSTANCES: [(u64, u64, u64, u64, u64); 2] = [
        // (modulus, order, base, logarithm, element)
        (2199023255867, 1099511627933, 9, 987654321, 1466426439827),
        (2199023258567, 1099511629283, 9, 1099511627775, 826684467692),
    ];

    #[test]
    #[cfg(feature = "parallel")]
    fn parallel() {
        // The walks are not reproducible with threads: the logarithm is what is checked, not how
        // it was reached.
        for (n, order, b, x, a) in PRIME_ORDER_INSTANCES {
            for threads in [2, 4] {
                assert_eq!(
                    discrete_log_pollard_rho_parallel(
                        &Integer::from(n),
                        &Integer::from(a),
                        &Integer::from(b),
                        Some(&Integer::from(order)),
                        threads,
                    )
                    .unwrap(),
                    x,
                    "log of {a} in base {b} modulo {n} over {threads} threads"
                );
            }
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn parallel_without_a_logarithm() {
        // 7 is not a power of 31 modulo 11: every thread walks its whole budget and stops, where
        // an unbounded search would never return.
        for threads in [0, 1, 2, 4, 8] {
            assert_eq!(
                discrete_log_pollard_rho_parallel(&11.into(), &7.into(), &31.into(), None, threads),
                Err(Error::LogDoesNotExist),
                "{threads} threads"
            );
        }
        // The same over GMP, with a modulus above `2**64`: `3` is not in the subgroup of order
        // 1009 that `b` generates, being of order `n - 1`.
        let order = Integer::from(1009u32);
        let (n, b) = big_modulus_and_base(&order);
        for threads in [1, 2, 4] {
            assert_eq!(
                discrete_log_pollard_rho_parallel(&n, &3.into(), &b, Some(&order), threads),
                Err(Error::LogDoesNotExist),
                "{threads} threads"
            );
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn parallel_thread_count_is_clamped() {
        // Counts no machine can spawn: a hundred thousand threads was a panic out of the scope
        // (`WouldBlock`), and `usize::MAX` a capacity overflow on the vector of seeds.
        for threads in [MAX_THREADS, MAX_THREADS + 1, 100_000, usize::MAX] {
            assert_eq!(
                discrete_log_pollard_rho_parallel(&11.into(), &7.into(), &31.into(), None, threads),
                Err(Error::LogDoesNotExist),
                "{threads} threads"
            );
        }
        // And a search over the cap still solves: the walks the clamp leaves are a whole search.
        assert_eq!(
            discrete_log_pollard_rho_parallel(
                &Integer::from(2147483783u32),
                &Integer::from(240127149u32),
                &9.into(),
                Some(&Integer::from(1073741891u32)),
                usize::MAX,
            )
            .unwrap(),
            12345678
        );
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn parallel_without_threads() {
        // No thread asked for, and one thread asked for, both run the sequential walk.
        for threads in [0, 1] {
            assert_eq!(
                discrete_log_pollard_rho_parallel(
                    &24567899.into(),
                    &(Integer::from(3).pow(333)),
                    &3.into(),
                    None,
                    threads,
                )
                .unwrap(),
                333,
                "{threads} threads"
            );
            assert_eq!(
                discrete_log_pollard_rho_parallel(
                    &227.into(),
                    &(Integer::from(3).pow(7)),
                    &5.into(),
                    None,
                    threads,
                )
                .unwrap(),
                132,
                "{threads} threads"
            );
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn parallel_modulus_above_a_word() {
        // The GMP path of the parallel search: the group lives on `BigRing`, whose elements are
        // the keys of the shared table.
        let order = Integer::from(1009u32);
        let (n, b) = big_modulus_and_base(&order);
        for x in [0u32, 1, 225, 1008] {
            let a = b.clone().pow_mod(&Integer::from(x), &n).unwrap();
            for threads in [2, 4] {
                assert_eq!(
                    discrete_log_pollard_rho_parallel(&n, &a, &b, Some(&order), threads).unwrap(),
                    x,
                    "log of {a} in base {b} modulo {n} over {threads} threads"
                );
            }
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn distinguished_points() {
        // A quarter of the bits of the order, and never more than the cap.
        assert_eq!(distinguished_bits(&Integer::from(4)), 1);
        assert_eq!(distinguished_bits(&(Integer::from(1) << 40u32)), 10);
        assert_eq!(
            distinguished_bits(&(Integer::from(1) << 400u32)),
            MAX_DISTINGUISHED_BITS
        );

        // One point in `2**bits` is distinguished, and the same point is distinguished in whichever
        // ring the group is held.
        let group = WordRing::new(1_000_003).unwrap();
        let big = BigRing::new(&Integer::from(1_000_003)).unwrap();
        for distinguished_bits in [1, 4, 8] {
            let search_group = Search {
                group: &group,
                exponents: &group,
                a: group.one(),
                b: group.one(),
                order: &Integer::from(1_000_003),
                multipliers: Vec::new(),
                budget: 0,
                distinguished_bits,
                shared: Shared {
                    table: Mutex::new(HashMap::new()),
                    answer: Mutex::new(None),
                    stop: AtomicBool::new(false),
                },
            };
            let mut met = 0usize;
            for value in 0..1_000_003u64 {
                if search_group.is_distinguished(&group.from_u64(value)) {
                    met += 1;
                    assert!(
                        big.residue_word(&big.from_integer(&Integer::from(value)))
                            .wrapping_mul(DISTINGUISHED_MULTIPLIER)
                            >> (u64::BITS - distinguished_bits)
                            == 0
                    );
                }
            }
            let expected = 1_000_003usize >> distinguished_bits;
            assert!(
                met > expected / 2 && met < expected * 2,
                "{met} of a million points distinguished at {distinguished_bits} bits, \
                 about {expected} expected"
            );
        }
    }
}
