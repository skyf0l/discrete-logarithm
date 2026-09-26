use std::{
    collections::HashMap,
    sync::{Arc, Mutex, PoisonError},
};

use primal::Sieve;
use rug::{Integer, rand::RandState};

use crate::{
    Error, check_modulus,
    modular::{BigRing, ModRing, WordRing},
};

/// Largest relation matrix the elimination is allowed to build, in bytes.
///
/// The elimination is dense: one row per prime of the factor base, of one entry per prime plus the
/// exponent of the base the relation stands for, and an entry is a machine word whenever the order
/// fits in one. That matrix, and not the running time, is what first puts an instance out of reach,
/// and it grows fast with the modulus: 836 primes and 5.3 MiB of rows for a modulus of 72 bits
/// (about 34 seconds of work), 2146 primes and 35 MiB at 88 bits, 17859 primes and 2.4 GiB at 128
/// bits, 49 GiB at 160 bits and 1.3 PiB at 256 bits. Above about 360 bits even the smoothness bound
/// is past `u32::MAX`.
///
/// 32 MiB is where the line is drawn. It is six times the rows of the largest instance ever
/// measured here, so nothing that was solvable is refused, and it is two orders of magnitude below
/// what a 128 bit modulus asks for, which is the size at which the search would otherwise spin over
/// candidates for a relation matrix no machine can hold. An order above `2**64` holds each entry in
/// a GMP integer instead of a word, several times more memory than this counts, which is the other
/// reason not to draw the line at what a machine could just about survive.
const MAX_RELATION_BYTES: u64 = 32 << 20;

/// Most primes the factor base may hold: the largest base [`MAX_RELATION_BYTES`] pays the rows of,
/// a base of `primes` holding up to `primes` rows of `primes + 1` elements.
const MAX_FACTORBASE: usize = {
    let element = std::mem::size_of::<u64>() as u64;
    let mut primes: u64 = 0;
    while (primes + 1) * (primes + 2) * element <= MAX_RELATION_BYTES {
        primes += 1;
    }
    primes as usize
};

/// Largest smoothness bound accepted: the 2048th prime, so that at most [`MAX_FACTORBASE`] primes
/// are below it.
///
/// The bound is what is tested, and before the sieve rather than after it: a modulus of 256 bits
/// asks for every prime below 71 million, whose sieve alone is tens of megabytes and seconds of
/// work, and the answer is going to be a refusal either way.
const MAX_BOUND: u32 = 17_863;

/// Number of partial relations held at once.
///
/// A partial relation is worth keeping only until another candidate meets the same large prime,
/// and how often that happens grows with the square of how many are held: a few thousand is where
/// the matches start on the instances this algorithm is chosen for. Each one costs the handful of
/// primes that divide its candidate, so the whole store stays well under a megabyte.
const MAX_PARTIALS: usize = 1 << 13;

/// Number of candidates one walk draws a starting point for.
///
/// A candidate drawn on its own is a whole powering, where the next one of a walk is a single
/// multiplication: a walk spreads the two powerings it starts from over this many candidates,
/// about an eighth of one each. A longer walk saves little more of that and lets a single starting
/// point decide most of the search instead, which the handful of relations a small factor base
/// needs leaves no room to recover from: an unlucky walk is dropped for a fresh one long before
/// the search gives up.
const WALK_STEPS: u32 = 16;

/// A residue about to be tested for smoothness.
enum Candidate {
    /// One that fits in a word, which is every residue of a modulus below `2**64`.
    Word(u64),
    /// One only GMP can hold.
    Big(Integer),
}

impl Candidate {
    /// The candidate `n`, taken as a word whenever it is one.
    ///
    /// Only the smoothness test exposed for the benchmarks and the tests of this module build one
    /// from an [`Integer`]: the algorithm itself gets its candidates from the ring the group
    /// arithmetic runs in, which knows whether they fit in a word.
    #[cfg(any(test, feature = "bench"))]
    fn of(n: Integer) -> Self {
        match n.to_u64() {
            Some(word) => Self::Word(word),
            None => Self::Big(n),
        }
    }
}

/// A ring the group arithmetic runs in, modulo `n`.
///
/// The one thing the index calculus needs of it beyond [`ModRing`] is the residue itself, to trial
/// divide: a ring holding it in a word hands it over without going through GMP.
trait GroupRing: ModRing {
    /// The residue of `x`, as the smoothness test divides it.
    fn candidate(&self, x: &Self::Elem) -> Candidate;
}

impl GroupRing for WordRing {
    fn candidate(&self, x: &u64) -> Candidate {
        Candidate::Word(self.to_u64(*x))
    }
}

impl GroupRing for BigRing {
    fn candidate(&self, x: &Integer) -> Candidate {
        match x.to_u64() {
            Some(word) => Candidate::Word(word),
            None => Candidate::Big(x.clone()),
        }
    }
}

/// A ring the relations are solved in, modulo the order.
///
/// The entries of a relation are exponents of a factorization and exponents of the base: small
/// numbers, which a ring holding a word-size order converts without building an [`Integer`].
// `from_word` takes a receiver because the modulus belongs to the ring, as `ModRing` documents.
#[allow(clippy::wrong_self_convention)]
trait OrderRing: ModRing {
    /// The element `x mod order`.
    fn from_word(&self, x: u64) -> Self::Elem;
}

impl OrderRing for WordRing {
    fn from_word(&self, x: u64) -> u64 {
        self.from_u64(x)
    }
}

impl OrderRing for BigRing {
    fn from_word(&self, x: u64) -> Integer {
        self.from_integer(&Integer::from(x))
    }
}

/// How a candidate factors over the factor base.
enum Smoothness {
    /// It factors completely: a relation, the exponents written to the buffer of the test.
    Smooth,
    /// It factors except for one prime cofactor, which is the value here: a partial relation, of
    /// use once a second candidate meets the same prime.
    LargePrime(u64),
    /// It does not factor, and nothing can be made of it.
    Rough,
}

/// The primes a candidate is divided by, and the bound that tells a large prime from a cofactor
/// nothing can be done with.
struct FactorBase<'p> {
    /// Every prime below the bound `B`, in increasing order.
    ///
    /// That all of them are there is what makes a cofactor below `B**2` a prime: such a cofactor
    /// has no factor below `B` left, and two factors at or above `B` would reach `B**2`.
    primes: &'p [u32],
    /// `B**2`, the largest cofactor that is still a single prime.
    large_prime_bound: u64,
}

impl<'p> FactorBase<'p> {
    /// The base of the primes below `bound`, which `primes` must hold all of.
    fn new(primes: &'p [u32], bound: u32) -> Self {
        Self {
            primes,
            large_prime_bound: u64::from(bound) * u64::from(bound),
        }
    }

    /// The base of exactly `primes`, the bound taken as small as they allow.
    ///
    /// What [`is_smooth`] needs: a caller passing a base of its own says nothing of the primes it
    /// left out, so a cofactor is only ever a large prime when it is below the square of the
    /// largest prime given.
    #[cfg(any(test, feature = "bench"))]
    fn of_primes(primes: &'p [u32]) -> Self {
        let bound = primes.last().map_or(0, |&last| last);
        Self::new(primes, bound)
    }

    /// The outcome of a candidate divided by every prime of the base, `left` being what the
    /// divisions left of it.
    fn verdict(&self, left: u64) -> Smoothness {
        if left == 1 {
            Smoothness::Smooth
        } else if left < self.large_prime_bound {
            // No prime of the base divides what is left and it is below `B**2`, so it is a single
            // prime: a second candidate meeting it cancels it.
            Smoothness::LargePrime(left)
        } else {
            Smoothness::Rough
        }
    }
}

/// The one factor base kept between calls, with the bound it was sieved for.
type FactorbaseCache = Mutex<Option<(u32, Arc<[u32]>)>>;

/// The factor base of the last bound sieved, held for the next call that asks for the same one.
///
/// Pohlig-Hellman solves one logarithm per prime power of the order, and every one of those calls
/// is modulo the same `n`, so every one of them sieves the same bound: about 5100 instructions for
/// the bound of 159 of a 28-bit modulus, and more the larger the modulus is. One base is all that
/// takes, which is also what keeps the memory bounded whatever a caller does.
static FACTORBASE: FactorbaseCache = Mutex::new(None);

/// The primes below `bound`, sieved unless [`FACTORBASE`] already holds exactly them.
fn factorbase(bound: u32) -> Arc<[u32]> {
    factorbase_from(&FACTORBASE, bound)
}

/// [`factorbase`] with the cache to look in, which the tests give one of their own.
///
/// A sieve of exactly that bound: the endless iterator of the primes sieves a fixed first segment
/// of its own, which costs more than the whole rest of a small instance.
fn factorbase_from(cache: &FactorbaseCache, bound: u32) -> Arc<[u32]> {
    // A panic while the lock is held leaves the cache as it was, which is still a base for the
    // bound stored with it: nothing about it is worth propagating a poisoning for.
    if let Some((cached, primes)) = &*cache.lock().unwrap_or_else(PoisonError::into_inner) {
        // The bound the primes were sieved for is stored with them, so a base is never handed out
        // for another one: a shorter base would make a smooth candidate look rough, and a longer
        // one would change which cofactor counts as a large prime.
        if *cached == bound {
            return Arc::clone(primes);
        }
    }

    // Sieved with the lock released: a call for another bound has nothing to wait for.
    let primes: Arc<[u32]> = Sieve::new(bound as usize)
        .primes_from(0)
        .take_while(|&p| p < bound as usize)
        .map(|p| p as u32)
        .collect();
    *cache.lock().unwrap_or_else(PoisonError::into_inner) = Some((bound, Arc::clone(&primes)));
    primes
}

/// Check if a number can be factored using the given factor base.
/// Returns the exponents vector if smooth, None otherwise.
///
/// The factor base must be primes in increasing order.
///
/// Nothing inside the crate calls it: the algorithm tests its own candidates against the base it
/// sieved. It is there for the benchmarks of the smoothness test, hence the feature gate.
#[cfg(any(test, feature = "bench"))]
pub fn is_smooth(n: Integer, factorbase: &[u32]) -> Option<Vec<u32>> {
    let mut exponents = vec![0u32; factorbase.len()];
    match smoothness(
        Candidate::of(n),
        &FactorBase::of_primes(factorbase),
        &mut exponents,
    ) {
        Smoothness::Smooth => Some(exponents),
        Smoothness::LargePrime(_) | Smoothness::Rough => None,
    }
}

/// How `candidate` factors over `base`, the exponent of each prime of the base written to the
/// matching entry of `exponents`.
fn smoothness(candidate: Candidate, base: &FactorBase<'_>, exponents: &mut [u32]) -> Smoothness {
    debug_assert_eq!(exponents.len(), base.primes.len(), "one exponent per prime");
    debug_assert!(
        base.primes.windows(2).all(|pair| pair[0] < pair[1]),
        "a factor base that is not increasing"
    );
    exponents.fill(0);
    match candidate {
        Candidate::Word(word) => divide_word(word, base, exponents, 0),
        Candidate::Big(big) => divide_big(big, base, exponents),
    }
}

/// [`smoothness`] of a candidate that fits in a word, the primes before `from` already divided
/// out of it and counted in `exponents`.
fn divide_word(
    mut left: u64,
    base: &FactorBase<'_>,
    exponents: &mut [u32],
    from: usize,
) -> Smoothness {
    // Zero is divisible by every prime, so dividing it out never ends. No power of the base is
    // zero unless the modulus makes it nilpotent, and such a power is not smooth either.
    if left == 0 {
        return Smoothness::Rough;
    }

    for (&prime, exponent) in base.primes[from..].iter().zip(&mut exponents[from..]) {
        // Nothing is left to divide.
        if left == 1 {
            break;
        }
        let prime = u64::from(prime);
        while left % prime == 0 {
            *exponent += 1;
            left /= prime;
        }
    }

    base.verdict(left)
}

/// [`smoothness`] of a candidate no word can hold.
fn divide_big(mut left: Integer, base: &FactorBase<'_>, exponents: &mut [u32]) -> Smoothness {
    if left == 0 {
        return Smoothness::Rough;
    }

    for index in 0..base.primes.len() {
        // A candidate that has come down to a word is divided by the word path from here: single
        // instructions per prime instead of a call into GMP for each of them.
        if let Some(word) = left.to_u64() {
            return divide_word(word, base, exponents, index);
        }
        let prime = base.primes[index];
        while left.is_divisible_u(prime) {
            exponents[index] += 1;
            left.div_exact_u_mut(prime);
        }
    }

    match left.to_u64() {
        Some(word) => base.verdict(word),
        // Above a word, so above `B**2` as well: not a large prime either.
        None => Smoothness::Rough,
    }
}

/// A candidate that factored over the factor base except for one large prime.
struct Partial<E> {
    /// The primes of the base that divide the candidate, each with its exponent.
    ///
    /// A candidate has a handful of them where a row of the elimination is dense, and these are
    /// held by the thousand: only the ones that are there are worth the memory.
    factors: Vec<(u32, u32)>,
    /// The exponent `x` the candidate is `b**x` of, modulo the order.
    x: E,
}

impl<E> Partial<E> {
    /// The partial relation of the candidate `b**x` whose factorization is `exponents`.
    fn new(exponents: &[u32], x: E) -> Self {
        Self {
            factors: exponents
                .iter()
                .enumerate()
                .filter(|&(_, &exponent)| exponent != 0)
                .map(|(index, &exponent)| (index as u32, exponent))
                .collect(),
            x,
        }
    }
}

/// The row of the elimination of a candidate: the exponent of each prime of the base, and the
/// exponent the candidate is a power of the base to.
///
/// A row `(e, x)` stands for the equation `sum(e_i * log(p_i)) = x`, which a relation satisfies
/// as it is and the relation of `a` satisfies with `log(a)` added to its right hand side.
fn dense_row<O: OrderRing>(exps: &O, exponents: &[u32], x: &O::Elem) -> Vec<O::Elem> {
    let mut row: Vec<O::Elem> = exponents
        .iter()
        .map(|&exponent| exps.from_word(u64::from(exponent)))
        .collect();
    row.push(x.clone());
    row
}

/// The relation of two candidates sharing the same large prime.
///
/// That prime divides each of them exactly once (it is below the square of the bound, so it is
/// prime, and its square would be above that bound), so the difference of the two rows is a
/// relation over the factor base alone: `b**(x - x_partial)` is the ratio of the two smooth parts.
fn combine<O: OrderRing>(
    exps: &O,
    partial: &Partial<O::Elem>,
    exponents: &[u32],
    x: &O::Elem,
) -> Vec<O::Elem> {
    let mut row = dense_row(exps, exponents, x);
    let last = row.len() - 1;
    for &(index, exponent) in &partial.factors {
        let term = exps.from_word(u64::from(exponent));
        exps.sub_assign(&mut row[index as usize], &term);
    }
    exps.sub_assign(&mut row[last], &partial.x);
    row
}

/// Subtracts from `row` the multiple of `pivot_row` that clears column `pivot`.
///
/// A stored relation is zero before its own column and holds 1 there, so the columns before
/// `pivot` are already what they must be and the pivot column clears itself. `factor` and `term`
/// are scratch elements: a ring whose elements live on the heap reuses their buffers instead of
/// allocating one per entry.
fn eliminate<O: ModRing>(
    exps: &O,
    row: &mut [O::Elem],
    pivot_row: &[O::Elem],
    pivot: usize,
    zero: &O::Elem,
    factor: &mut O::Elem,
    term: &mut O::Elem,
) {
    factor.clone_from(&row[pivot]);
    for (entry, coefficient) in row.iter_mut().zip(pivot_row).skip(pivot) {
        // A relation in reduced form is zero in every column another pivot holds, and that is most
        // of them once the elimination is under way.
        if coefficient == zero {
            continue;
        }
        term.clone_from(coefficient);
        exps.mul_assign(term, factor);
        exps.sub_assign(entry, term);
    }
}

/// Index Calculus algorithm for computing the discrete logarithm of `a` in base `b` modulo `n`.
///
/// The group order must be given and prime. It is not suitable for small orders
/// and the algorithm might fail to find a solution in such situations.
///
/// This algorithm is particularly efficient for large prime orders when
/// exp(2*sqrt(log(n)*log(log(n)))) < sqrt(order).
///
/// A modulus that is not positive is refused with [`Error::InvalidModulus`]. Modulo 1 every residue
/// is 0, so the logarithm is 0. A modulus whose factor base needs more than 32 MiB of relation
/// matrix is refused with [`Error::OutOfReach`], which is from about 88 bits on: the matrix of a 128
/// bit modulus is 2.4 GiB, and a search that cannot hold its own matrix would only spin over
/// candidates for as long as it was left to.
///
/// [`Error::LogDoesNotExist`] from here means "no logarithm found within the budget of the search",
/// not a proof that none exists: the relations are looked for a bounded number of times, and an
/// instance the search gave up on comes back as that error. It is also what a base or a target that
/// is not invertible modulo `n` gets, the walk looking for the first relation reaching zero, which
/// no multiplication leaves. The logarithm is verified before it is returned, so a result is never
/// wrong; only a failure can be.
///
/// # Examples
///
/// ```
/// use discrete_logarithm::discrete_log_index_calculus;
/// use rug::Integer;
///
/// let n = Integer::from(24570203447_u64);
/// let a = Integer::from(23859756228_u64);
/// let b = Integer::from(2);
/// let order = Integer::from(12285101723_u64);
///
/// let x = discrete_log_index_calculus(&n, &a, &b, Some(&order)).unwrap();
/// assert_eq!(x, Integer::from(4519867240_u64));
/// ```
///
/// If the order of the group is known, it must be passed as `order`.
pub fn discrete_log_index_calculus(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    solve(n, a, b, order, &mut RandState::new())
}

/// Index Calculus algorithm for computing the discrete logarithm of `a` in base `b` modulo `n`.
///
/// Same as [`discrete_log_index_calculus`], with the random generator seeded with `seed`: the
/// same seed always looks for the relations in the same order, and a retry with another seed
/// looks for others.
///
/// The preconditions and the meaning of every error are the ones of
/// [`discrete_log_index_calculus`].
pub fn discrete_log_index_calculus_with_seed(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    seed: u64,
) -> Result<Integer, Error> {
    let mut rand_state = RandState::new();
    rand_state.seed(&Integer::from(seed));
    solve(n, a, b, order, &mut rand_state)
}

/// The smoothness bound `B` of the factor base of `n`, the heuristic of the sympy implementation:
/// `B = exp(0.5 * sqrt(log(n) * log(log(n))) * (1 + 1/log(log(n))))`.
///
/// Not a number for a modulus below 3, whose double logarithm is not one: every comparison with
/// `NaN` is false, so such a modulus is neither refused as too large nor accepted, and the empty
/// factor base it leads to is what turns it away.
fn smoothness_bound(n: &Integer) -> f64 {
    let log_n = n.to_f64().ln();
    let log_log_n = log_n.ln();
    (0.5 * (log_n * log_log_n).sqrt() * (1.0 + 1.0 / log_log_n)).exp()
}

fn solve(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error> {
    // The one entry point that used to reach its arithmetic without looking at the modulus: it came
    // out right only through a chain of coincidences (the logarithm of a negative number is `NaN`,
    // `NaN as u32` is 0, a bound of 0 sieves an empty factor base, an empty base is refused).
    check_modulus(n)?;
    // Modulo 1 every residue is 0, so every exponent is a logarithm and 0 is the smallest one.
    if *n == 1 {
        return Ok(Integer::new());
    }

    let order = match order {
        Some(order) => order,
        None => return Err(Error::LogDoesNotExist),
    };

    let b_bound = smoothness_bound(n);
    // A factor base whose relation matrix is past the budget: out of reach, which says nothing about
    // the logarithm. A modulus below 3, whose bound is not even a number, falls through instead: the
    // comparison is false for `NaN`, and the empty factor base below is what refuses it.
    if b_bound > f64::from(MAX_BOUND) {
        return Err(Error::OutOfReach);
    }
    let b_bound = b_bound as u32;

    // Compute the factorbase - all primes up to B (exclusive, matching sympy's primerange(B))
    let factorbase = factorbase(b_bound);
    // What the bound above is there to guarantee, so the relation matrix stays within
    // `MAX_RELATION_BYTES`: cheap enough to check that the two constants agree.
    debug_assert!(
        factorbase.len() <= MAX_FACTORBASE,
        "a factor base of {} primes, past the {MAX_FACTORBASE} the budget pays for",
        factorbase.len()
    );

    // No base at all, which is also how every modulus too small to have one is refused: the bound
    // of such a modulus is not even a number.
    if factorbase.is_empty() {
        return Err(Error::LogDoesNotExist);
    }
    let base = FactorBase::new(&factorbase, b_bound);

    // Maximum number of tries to find a relation
    let max_tries = 5 * u128::from(b_bound) * u128::from(b_bound);

    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);

    // The only exponent a group of this order has is 0, and no exponent can be drawn below it
    // either. A non positive order has no exponent at all.
    if *order < 2 {
        return if *order == 1 && a == 1 {
            Ok(Integer::new())
        } else {
            Err(Error::LogDoesNotExist)
        };
    }

    // The word-size ring whenever the modulus fits in one, GMP otherwise.
    match WordRing::from_modulus(n) {
        Some(group) => with_group(&group, &a, &b, order, &base, max_tries, rand_state),
        None => {
            let group = BigRing::new(n).expect("a modulus a factor base was built for");
            with_group(&group, &a, &b, order, &base, max_tries, rand_state)
        }
    }
}

/// [`solve`] with the group ring chosen, the ring of the relations left to choose.
fn with_group<G: GroupRing>(
    group: &G,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    base: &FactorBase<'_>,
    max_tries: u128,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error> {
    // The relation arithmetic is modulo the order, which is its own ring: an order of word size
    // makes a row of the elimination a row of words, whatever the modulus is.
    match WordRing::from_modulus(order) {
        Some(exps) => run(group, &exps, a, b, order, base, max_tries, rand_state),
        None => {
            let exps = BigRing::new(order).expect("an order above 1");
            run(group, &exps, a, b, order, base, max_tries, rand_state)
        }
    }
}

/// [`solve`] in `group` for the group arithmetic and `exps` for the relations.
#[allow(clippy::too_many_arguments)]
fn run<G: GroupRing, O: OrderRing>(
    group: &G,
    exps: &O,
    a: &Integer,
    b: &Integer,
    order: &Integer,
    base: &FactorBase<'_>,
    max_tries: u128,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error> {
    let lf = base.primes.len();
    let zero = exps.zero();
    let a = group.from_integer(a);
    let b = group.from_integer(b);
    let one = group.one();
    // The exponents of one candidate, reused by every smoothness test instead of allocated per
    // candidate: the tests are counted in thousands per relation.
    let mut exponents = vec![0u32; lf];

    // First, find a relation for a: an `x` with `a * b**x` smooth over the factorbase.
    let mut relationa = None;
    let mut abx = a.clone();
    let group_zero = group.zero();
    for x in 0..order.to_u64().unwrap_or(u64::MAX) {
        // Nothing multiplies a zero residue into anything else, and zero is not smooth: every
        // candidate left is this one, and walking the rest of the order only takes longer to say so.
        // A target of 0 starts here, and a base that is not invertible modulo `n` arrives here.
        if abx == group_zero {
            break;
        }

        if abx == one {
            return Ok((order.clone() - x) % order);
        }

        if let Smoothness::Smooth = smoothness(group.candidate(&abx), base, &mut exponents) {
            relationa = Some(dense_row(exps, &exponents, &exps.from_word(x)));
            break;
        }

        // A base of 1 leaves the candidate where it is, so every step of the walk would test the
        // residue just tested: one test is all there is to do. Nothing sends a base of 1 here on
        // purpose, and an order of `10**12` made the walk minutes long.
        if b == one {
            break;
        }

        group.mul_assign(&mut abx, &b);
    }

    let Some(mut relationa) = relationa else {
        return Err(Error::LogDoesNotExist);
    };

    // Now find relations for the factorbase elements: one pivot per prime, and the relation of `a`
    // reduced by each new one until no unknown is left in it.
    let mut relations: Vec<Option<Vec<O::Elem>>> = vec![None; lf];
    // Candidates that left one large prime, indexed by it: two of them make one relation.
    let mut partials: HashMap<u64, Partial<O::Elem>> = HashMap::new();
    // Scratch elements of the elimination, allocated once for the whole run.
    let mut factor = exps.zero();
    let mut term = exps.zero();

    // Pivots the elimination holds, one per relation that added information. Sympy counts every
    // relation it finds against a budget of `3 * lf` of them, so a relation already implied by the
    // ones before it, which the elimination throws away, still consumes the budget: that is what
    // makes a small instance give up while candidates are still coming. Only a relation that fills a
    // pivot brings the answer closer, and `lf` of them are all a target can need, the relation of
    // `a` being solved as soon as the pivots cover the primes left in it.
    let mut rank = 0;
    // Candidates drawn since the last relation that added information: what the search gives up on,
    // and, the rank being bounded by `lf`, the only limit that ever stops it.
    let mut failures = 0u128;
    let order_minus_1 = Integer::from(order - 1u32);

    // The walk the candidates come from: `bx` is `b**x`, and the next candidate is one
    // multiplication by `stride` away, where a powering of its own would be most of what a
    // candidate costs. The exponent of the stride is added to `x` at the same time, so the exponent
    // of every candidate is known without ever solving for it.
    let mut bx = group.zero();
    let mut stride = group.one();
    let mut x = exps.zero();
    let mut stride_exponent = exps.zero();
    // Nothing is walked yet, so the first candidate draws a starting point.
    let mut step = WALK_STEPS;

    while rank < lf && failures < max_tries {
        if step < WALK_STEPS {
            group.mul_assign(&mut bx, &stride);
            exps.add_assign(&mut x, &stride_exponent);
            step += 1;
        } else {
            // A fresh starting point and a fresh stride, both exponents in `[1, order - 1]`, and
            // the only two powerings of the walk: a square and multiply over words for a modulus
            // below `2**64`, GMP's own powering above it.
            //
            // The stride is drawn and not taken as `b` itself, although that would save one of the
            // two: the candidates of a walk of step `b` are `b` times one another, so when `b` is
            // smooth (base 2, with 2 in the factor base) all the relations of a walk follow from
            // the first one and the single equation `log(b) = 1`, and the elimination never fills
            // its pivots.
            let start = order_minus_1.clone().random_below(rand_state) + 1u32;
            let stride_start = order_minus_1.clone().random_below(rand_state) + 1u32;
            bx = group.pow(&b, &start);
            stride = group.pow(&b, &stride_start);
            x = exps.from_integer(&start);
            stride_exponent = exps.from_integer(&stride_start);
            step = 1;
        }

        // Try to factor it over the factorbase
        let mut relation = match smoothness(group.candidate(&bx), base, &mut exponents) {
            Smoothness::Smooth => dense_row(exps, &exponents, &x),
            Smoothness::LargePrime(prime) => match partials.get(&prime) {
                // The large prime cancels between the two, which leaves a relation over the
                // factor base: several times more of them for the same candidates.
                Some(partial) => combine(exps, partial, &exponents, &x),
                None => {
                    // Past the limit the store is left as it is: the partials already there
                    // are as likely to be matched as any new one.
                    if partials.len() < MAX_PARTIALS {
                        partials.insert(prime, Partial::new(&exponents, x.clone()));
                    }
                    failures += 1;
                    continue;
                }
            },
            Smoothness::Rough => {
                failures += 1;
                continue;
            }
        };

        // Gaussian elimination step: clear every column that already has a pivot, and keep the
        // first column that has none.
        let mut index = lf;
        for (column, pivot) in relations.iter().enumerate() {
            if relation[column] == zero {
                continue;
            }
            match pivot {
                Some(pivot) => eliminate(
                    exps,
                    &mut relation,
                    pivot,
                    column,
                    &zero,
                    &mut factor,
                    &mut term,
                ),
                None => {
                    if index == lf {
                        index = column;
                    }
                }
            }
        }

        // No new information: every prime of the relation already has a pivot, so the elimination
        // has nothing to keep. It counts as a candidate that brought nothing, which is what makes
        // the budget below one of informative relations.
        if index == lf {
            failures += 1;
            continue;
        }

        // Normalize the relation, so that eliminating with it is a multiplication by the entry to
        // clear.
        let Some(inverse) = exps.invert(&relation[index]) else {
            // Only an order that is not prime, which this algorithm is not for, has an entry that
            // cannot be inverted: drop the relation instead of failing on it.
            failures += 1;
            continue;
        };
        for entry in relation.iter_mut().skip(index) {
            exps.mul_assign(entry, &inverse);
        }
        relations[index] = Some(relation);
        rank += 1;
        failures = 0;

        // Reduce relationa with the new relation, and with the ones it opens the way to.
        for (column, pivot) in relations.iter().enumerate() {
            if relationa[column] == zero {
                continue;
            }
            match pivot {
                // A nonzero entry with no pivot: nothing more can be eliminated.
                None => break,
                Some(pivot) => eliminate(
                    exps,
                    &mut relationa,
                    pivot,
                    column,
                    &zero,
                    &mut factor,
                    &mut term,
                ),
            }
        }

        // Check if all unknowns are eliminated
        if relationa[..lf].iter().all(|entry| *entry == zero) {
            // What is left of the relation of `a` is `log(a) = -relationa[lf]`.
            let x = exps.to_integer(&exps.sub(&zero, &relationa[lf]));

            // Verify the result: the relations only hold if the order given really is the order of
            // the base, and a wrong one gives a wrong logarithm.
            if group.pow(&b, &x) == a {
                return Ok(x);
            }

            return Err(Error::LogDoesNotExist);
        }
    }

    Err(Error::LogDoesNotExist)
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
    fn index_calculus() {
        // Test case from sympy documentation
        assert_eq!(
            discrete_log_index_calculus(
                &int("24570203447"),
                &int("23859756228"),
                &2.into(),
                Some(&int("12285101723"))
            )
            .unwrap(),
            int("4519867240")
        );
    }

    #[test]
    fn index_calculus_small() {
        // Small test cases
        assert_eq!(
            discrete_log_index_calculus(
                &587.into(),
                &(Integer::from(2).pow(9)),
                &2.into(),
                Some(&293.into())
            )
            .unwrap(),
            9
        );
    }

    #[test]
    fn sympy_cases() {
        for (n, a, b, order, x) in [
            ("983", "948", "2", "491", "183"),
            ("633383", "21794", "2", "316691", "68048"),
            ("941762639", "68822582", "2", "470881319", "338029275"),
            (
                "999231337607",
                "888188918786",
                "2",
                "499615668803",
                "142811376514",
            ),
            // A base much larger than the primes of the factorbase.
            (
                "47747730623",
                "19410045286",
                "43425105668",
                "645239603",
                "590504662",
            ),
        ] {
            assert_eq!(
                discrete_log_index_calculus(&int(n), &int(a), &int(b), Some(&int(order))).unwrap(),
                int(x),
                "n = {n}"
            );
        }
    }

    #[test]
    fn seeded() {
        // The same seed always looks for the same relations, and the result does not depend on it.
        for seed in [0, 1, 42] {
            assert_eq!(
                discrete_log_index_calculus_with_seed(
                    &int("941762639"),
                    &int("68822582"),
                    &2.into(),
                    Some(&int("470881319")),
                    seed
                )
                .unwrap(),
                int("338029275")
            );
        }
    }

    #[test]
    fn every_walk_solves() {
        // Every starting point leads to the same logarithm, base 2 included: the candidates of a
        // walk must not be multiples of one another, which taking the base itself as the stride
        // would make them whenever the base is smooth. Every relation of such a walk follows from
        // its first one and the single equation `log(b) = 1`, and the relation budget goes on them.
        for seed in 0..24 {
            assert_eq!(
                discrete_log_index_calculus_with_seed(
                    &int("941762639"),
                    &int("68822582"),
                    &2.into(),
                    Some(&int("470881319")),
                    seed
                )
                .unwrap(),
                int("338029275"),
                "seed = {seed}"
            );
        }
    }

    #[test]
    fn factorbase_cache() {
        // A cache of its own, so that the other tests running beside this one cannot evict what it
        // is about to read back.
        let cache = FactorbaseCache::new(None);

        let first = factorbase_from(&cache, 30);
        assert_eq!(*first, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]);
        // The same bound comes back as the very same base, not sieved again.
        let second = factorbase_from(&cache, 30);
        assert!(Arc::ptr_eq(&first, &second));

        // Another bound is sieved, and never answered with the base of the previous one.
        let shorter = factorbase_from(&cache, 10);
        assert_eq!(*shorter, [2, 3, 5, 7]);
        // And back to the first bound, which the one base held cannot answer any more.
        let again = factorbase_from(&cache, 30);
        assert_eq!(*again, *first);
        assert!(!Arc::ptr_eq(&first, &again));

        // The bounds no prime is below, which is how a modulus too small is refused.
        assert!(factorbase_from(&cache, 2).is_empty());
        assert!(factorbase_from(&cache, 0).is_empty());
    }

    #[test]
    fn no_order() {
        // The order is required.
        assert_eq!(
            discrete_log_index_calculus(&int("24570203447"), &int("23859756228"), &2.into(), None),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn smoothness() {
        let factorbase = [2, 3, 5, 7];
        assert_eq!(
            is_smooth(Integer::from(2 * 2 * 3 * 7 * 7), &factorbase),
            Some(vec![2, 1, 0, 2])
        );
        assert_eq!(is_smooth(Integer::from(1), &factorbase), Some(vec![0; 4]));
        assert_eq!(is_smooth(Integer::from(11), &factorbase), None);
        assert_eq!(is_smooth(Integer::from(2 * 11), &factorbase), None);
    }

    #[test]
    fn smoothness_of_zero() {
        // Zero is divisible by every prime: dividing it out never ends, where sympy loops forever.
        let factorbase = [2, 3, 5, 7];
        assert_eq!(is_smooth(Integer::new(), &factorbase), None);
        // A candidate GMP holds, of the same kind: `-2**80` divides down to `-1`.
        assert_eq!(is_smooth(Integer::from(-1) << 80u32, &factorbase), None);
        // And the way it is reached from the outside: `a = 0` is no power of the base.
        assert_eq!(
            discrete_log_index_calculus(&int("983"), &Integer::new(), &2.into(), Some(&int("491"))),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn large_prime_cofactors() {
        // Primes below 30, and candidates whose cofactor is a single prime above them: what the
        // partial relations are made of.
        let factorbase: Vec<u32> = Sieve::new(30)
            .primes_from(0)
            .take_while(|&p| p < 30)
            .map(|p| p as u32)
            .collect();
        let base = FactorBase::new(&factorbase, 30);
        let mut exponents = vec![0u32; factorbase.len()];

        let smooth = 2 * 3 * 23;
        assert!(matches!(
            super::smoothness(Candidate::of(Integer::from(smooth)), &base, &mut exponents),
            Smoothness::Smooth
        ));
        assert_eq!(exponents[0], 1);
        assert_eq!(exponents[1], 1);
        // 23 is the prime of index 8.
        assert_eq!(exponents[8], 1);

        // A prime above the base but below the square of its bound: a partial relation.
        assert!(matches!(
            super::smoothness(Candidate::of(Integer::from(4 * 31)), &base, &mut exponents),
            Smoothness::LargePrime(31)
        ));
        assert_eq!(exponents[0], 2);

        // Two primes above the base, which is what the bound tells apart from one: nothing can be
        // done with it.
        assert!(matches!(
            super::smoothness(Candidate::of(Integer::from(31 * 37)), &base, &mut exponents),
            Smoothness::Rough
        ));
        // Neither can with a single prime at or above that square.
        assert!(matches!(
            super::smoothness(Candidate::of(Integer::from(907)), &base, &mut exponents),
            Smoothness::Rough
        ));

        // A candidate too large for a word, which comes down to one as it is divided.
        let big = Integer::from(2).pow(80u32) * 31u32;
        assert!(matches!(
            super::smoothness(Candidate::of(big), &base, &mut exponents),
            Smoothness::LargePrime(31)
        ));
        assert_eq!(exponents[0], 80);
    }

    #[test]
    fn partial_relations_combine() {
        // Two candidates sharing a large prime: their difference is a relation over the base
        // alone, the large prime cancelled.
        let exps = WordRing::new(1000003).unwrap();
        let first = Partial::new(&[1, 0, 2], exps.from_word(7));
        let second = [0u32, 3, 1];
        let row = combine(&exps, &first, &second, &exps.from_word(20));
        let residues: Vec<Integer> = row.iter().map(|entry| exps.to_integer(entry)).collect();
        assert_eq!(
            residues,
            vec![
                Integer::from(1000003 - 1),
                Integer::from(3),
                Integer::from(1000003 - 1),
                Integer::from(13)
            ]
        );
    }

    #[test]
    fn new_instances() {
        // Instances of the size this algorithm is used on, every logarithm verified against
        // sympy's `_discrete_log_index_calculus`: `n = 2 * order + 1` is a safe prime, and the
        // base is a square, of order the large prime factor of `n - 1`.
        for (n, a, b, order, x) in [
            // 35, 38 and 41 bits, the order one bit below each.
            (
                "21562317419",
                "12917980622",
                "4",
                "10781158709",
                "10443968216",
            ),
            (
                "168477881123",
                "31620352101",
                "4",
                "84238940561",
                "66857962857",
            ),
            (
                "1188783683903",
                "850486924030",
                "4",
                "594391841951",
                "571046147851",
            ),
            // An order that is a word-size prime just above `2**40`.
            (
                "2843340667127",
                "1845576079565",
                "4",
                "1421670333563",
                "20782261317",
            ),
            // A base far above every prime of the factorbase.
            (
                "1188783683903",
                "697702363150",
                "130758869943",
                "594391841951",
                "503471300153",
            ),
        ] {
            assert_eq!(
                discrete_log_index_calculus(&int(n), &int(a), &int(b), Some(&int(order))).unwrap(),
                int(x),
                "n = {n}"
            );
            // The same logarithm whatever the relations found are.
            for seed in [0, 7, 0xD15C] {
                assert_eq!(
                    discrete_log_index_calculus_with_seed(
                        &int(n),
                        &int(a),
                        &int(b),
                        Some(&int(order)),
                        seed
                    )
                    .unwrap(),
                    int(x),
                    "n = {n}, seed = {seed}"
                );
            }
        }
    }

    #[test]
    fn degenerate_orders() {
        // Orders no exponent can be drawn below, where sympy panics on an empty range.
        for order in [-1, 0, 2] {
            assert_eq!(
                discrete_log_index_calculus(
                    &int("941762639"),
                    &int("68822582"),
                    &2.into(),
                    Some(&order.into())
                ),
                Err(Error::LogDoesNotExist)
            );
        }
        // The only logarithm a group of order 1 has.
        assert_eq!(
            discrete_log_index_calculus(&int("941762639"), &1.into(), &2.into(), Some(&1.into()))
                .unwrap(),
            0
        );
    }

    #[test]
    fn invalid_moduli() {
        // Checked like every other entry point, instead of being left to the coincidences of a
        // logarithm of a negative number.
        for n in [-4, -1, 0] {
            assert_eq!(
                discrete_log_index_calculus(&n.into(), &1.into(), &2.into(), Some(&491.into())),
                Err(Error::InvalidModulus),
                "n = {n}"
            );
        }
        // Modulo 1 every residue is 0, so the logarithm is 0.
        assert_eq!(
            discrete_log_index_calculus(&1.into(), &1.into(), &2.into(), Some(&491.into()))
                .unwrap(),
            0
        );
        assert_eq!(
            discrete_log_index_calculus_with_seed(&1.into(), &7.into(), &5.into(), None, 3)
                .unwrap(),
            0
        );
        // A modulus with no factor base at all: 2 is too small for its bound to even be a number,
        // and the empty base is what turns it away.
        assert!(smoothness_bound(&2.into()).is_nan());
        assert_eq!(
            discrete_log_index_calculus(&2.into(), &1.into(), &2.into(), Some(&491.into())),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn walk_reaching_zero_stops() {
        // `b = 0` modulo 4: the walk looking for the relation of `a` reaches zero at its first step
        // and nothing multiplies it out of there. It used to run over the whole order to learn it,
        // which is minutes at `10**8` and hours at `10**12`.
        for order in [
            Integer::from(1_000_000),
            Integer::from(100_000_000),
            Integer::from(1_000_000_000_000u64),
            Integer::from(1) << 80u32,
        ] {
            assert_eq!(
                discrete_log_index_calculus(&4.into(), &3.into(), &0.into(), Some(&order)),
                Err(Error::LogDoesNotExist),
                "order {order}"
            );
        }
        // The same for a base that is a zero divisor rather than zero: `2**2 = 0` modulo 4.
        assert_eq!(
            discrete_log_index_calculus(
                &4.into(),
                &3.into(),
                &2.into(),
                Some(&Integer::from(1_000_000_000_000u64))
            ),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn walk_of_a_base_of_one_stops() {
        // `b = 1` modulo 4, as 1 and -7 both are: every candidate of the walk is `a` itself, so one
        // test settles it. The walk used to run over the whole order instead.
        for b in [1, -7, 5] {
            assert_eq!(
                discrete_log_index_calculus(
                    &4.into(),
                    &(-1).into(),
                    &b.into(),
                    Some(&Integer::from(1_000_000_000_000u64))
                ),
                Err(Error::LogDoesNotExist),
                "b = {b}"
            );
        }
        // And the one logarithm a base of 1 has.
        assert_eq!(
            discrete_log_index_calculus(
                &4.into(),
                &1.into(),
                &1.into(),
                Some(&Integer::from(1_000_000_000_000u64))
            )
            .unwrap(),
            0
        );
    }

    #[test]
    fn factor_base_cap() {
        // The cap is the largest factor base the relation matrix budget pays for.
        assert_eq!(factorbase(MAX_BOUND).len(), MAX_FACTORBASE);
        assert!(factorbase(MAX_BOUND + 1).len() > MAX_FACTORBASE);
    }

    #[test]
    fn modulus_out_of_reach() {
        // Safe primes of 72, 86 and 88 bits: the first two are accepted (836 primes and 5.3 MiB of
        // rows, 1916 and 28 MiB), the third is not (2146 and 35 MiB). Solving the accepted ones
        // takes half a minute and several minutes, which no test should wait for: the bound they are
        // accepted by is what is checked here, and the refusal of the third end to end below.
        for accepted in [
            int("4169048128772662588679"),
            int("72398380214093650144774079"),
        ] {
            assert!(
                smoothness_bound(&accepted) <= f64::from(MAX_BOUND),
                "{accepted} of {} bits refused",
                accepted.significant_bits()
            );
        }

        // The refusal is immediate, and says nothing about the logarithm: each of these exists, `4`
        // generating the subgroup of order `(n - 1) / 2` of a safe prime.
        for (n, a) in [
            (
                int("245827555352776502451081623"),
                int("191727569291595824293007544"),
            ),
            // 128 bits, whose 17859 primes would need 2.4 GiB of rows.
            (
                int("279305127720133152706028210119365912959"),
                int("269390257626015088577421517053028557938"),
            ),
            // And far above, where the smoothness bound itself is past `u32::MAX`: sieving it is
            // what the bound check keeps out.
            ((Integer::from(1) << 400u32).next_prime(), Integer::from(4)),
        ] {
            assert!(
                smoothness_bound(&n) > f64::from(MAX_BOUND),
                "{n} accepted by the bound"
            );
            let order = Integer::from(&n - 1u32) / 2u32;
            assert_eq!(
                discrete_log_index_calculus(&n, &a, &4.into(), Some(&order)),
                Err(Error::OutOfReach),
                "n of {} bits",
                n.significant_bits()
            );
        }
    }
}
