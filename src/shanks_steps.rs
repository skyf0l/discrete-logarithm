use std::hash::{Hash, Hasher};

use rug::Integer;

use crate::{
    Error, check_modulus,
    modular::{BigRing, ModRing, WordRing},
    n_order,
};

/// Largest order the table of baby steps is built for, a deliberately conservative cap.
///
/// What it costs and why it is where it is: see [`discrete_log_shanks_steps`].
pub const MAX_ORDER: u64 = 1_000_000_000_000u64;

/// Baby steps taken between two giant steps, `c` below.
///
/// `k` rounds of `c` baby steps cover about `c * k**2 / 2` exponents for `(c + 2) * k`
/// multiplications, so a logarithm `x` costs about `(c + 2) * sqrt(2 * x / c)` of them. Once the
/// baby steps stop, at the square root `m` of the order, the search is a plain baby-step
/// giant-step one, of which the rounds before have cost about `2.5 * m / c` multiplications more:
/// sixteen keeps a large logarithm within a tenth of what it costs without interleaving, and a
/// small one within a small multiple of `sqrt(x)`.
const BABY_STEPS_PER_ROUND: u64 = 16;

/// Exponent marking a free slot of [`BabySteps`].
///
/// No baby step reaches it: their exponents stop at the square root of [`MAX_ORDER`].
const FREE: u32 = u32::MAX;

/// Number of slots of a new [`BabySteps`]: what the first rounds of a small logarithm need,
/// and small enough that starting below it would save nothing.
const INITIAL_SLOTS: usize = 16;

/// Multiplier of the hash: the odd number closest to `2**64 / phi`, whose multiples spread over
/// the whole word, the high bits of the product being the mixed ones.
const HASH_MULTIPLIER: u64 = 0x9E37_79B9_7F4A_7C15;

/// Hasher of [`BabySteps`]: one multiplication per word of the key.
///
/// The keys are residues, already spread over the modulus, and nobody but this algorithm chooses
/// them: they only have to be mixed into the high bits of a word, where SipHash (what a `HashMap`
/// hashes with) runs a dozen rounds per key to also resist an adversary picking them.
#[derive(Default)]
struct MulShift(u64);

impl Hasher for MulShift {
    fn finish(&self) -> u64 {
        self.0
    }

    #[inline]
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0 ^ x).wrapping_mul(HASH_MULTIPLIER);
    }

    fn write(&mut self, bytes: &[u8]) {
        // The way a key made of anything else than words reaches the hash: GMP hashes a bignum
        // limb by limb, a `u64` in a single `write_u64`.
        for chunk in bytes.chunks(8) {
            let mut word = [0u8; 8];
            word[..chunk.len()].copy_from_slice(chunk);
            self.write_u64(u64::from_le_bytes(word));
        }
    }
}

/// The hash of `key`.
#[inline]
fn hash_key<E: Hash>(key: &E) -> u64 {
    let mut hasher = MulShift::default();
    key.hash(&mut hasher);
    hasher.finish()
}

/// The baby steps taken so far, each mapped to its exponent.
///
/// Open addressing with linear probing over a table of exponents, next to the baby steps
/// themselves in a dense vector: a free slot is four bytes to write, where a slot holding a key
/// would cost a copy of one (and a GMP allocation per slot for a modulus above `2**64`), and a
/// `HashMap<Integer, u64>` costs an allocation per entry and a SipHash per lookup. The table starts
/// with [`INITIAL_SLOTS`] slots and doubles, so a small logarithm never pays for a table sized for
/// the whole order.
struct BabySteps<E> {
    /// `b**0, b**1, ...`, the exponent of a baby step being its position.
    babies: Vec<E>,
    /// The exponent each slot holds, [`FREE`] marking a free one. Their number is a power of two,
    /// which is what masking a hash into an index needs.
    slots: Vec<u32>,
    /// `64 - log2(slots.len())`: a hash is mixed into its high bits, an index comes from there.
    shift: u32,
}

impl<E: Eq + Hash> BabySteps<E> {
    /// An empty table.
    fn new() -> Self {
        Self {
            babies: Vec::new(),
            slots: vec![FREE; INITIAL_SLOTS],
            shift: u64::BITS - INITIAL_SLOTS.trailing_zeros(),
        }
    }

    /// The exponent of `key`, `None` when it is not one of the baby steps.
    #[inline]
    fn get(&self, key: &E) -> Option<u32> {
        let slots = self.slots.as_slice();
        // An index masked with the number of slots less one is inside them, whatever the hash is:
        // the compiler sees it from the slice itself and leaves the bounds checks out.
        let mask = slots.len() - 1;
        let mut index = (hash_key(key) >> self.shift) as usize & mask;
        loop {
            let exponent = slots[index];
            if exponent == FREE {
                return None;
            }
            if &self.babies[exponent as usize] == key {
                return Some(exponent);
            }
            index = (index + 1) & mask;
        }
    }

    /// Adds `baby`, the baby step of the exponent that follows the last one.
    #[inline]
    fn push(&mut self, baby: E) {
        let exponent = self.babies.len() as u32;
        debug_assert!(exponent != FREE, "exponent of a free slot");
        self.babies.push(baby);
        self.index(exponent);
        // At most half of the slots are occupied: a fuller table makes linear probing walk over
        // long runs of occupied ones, and that walk is what a lookup and an insertion cost.
        if self.babies.len() * 2 >= self.slots.len() {
            self.grow();
        }
    }

    /// Puts `exponent` in the slots, unless its baby step is already there.
    ///
    /// Baby steps are added by increasing exponent, and a repeated power of the base is a power the
    /// smaller exponent already reaches: only the smallest solution is wanted.
    #[inline]
    fn index(&mut self, exponent: u32) {
        let shift = self.shift;
        let babies = self.babies.as_slice();
        let slots = self.slots.as_mut_slice();
        let mask = slots.len() - 1;
        let baby = &babies[exponent as usize];
        let mut index = (hash_key(baby) >> shift) as usize & mask;
        loop {
            let slot = slots[index];
            if slot == FREE {
                slots[index] = exponent;
                return;
            }
            if &babies[slot as usize] == baby {
                return;
            }
            index = (index + 1) & mask;
        }
    }

    /// Doubles the slots, the baby steps indexed into them again.
    fn grow(&mut self) {
        self.slots = vec![FREE; self.slots.len() * 2];
        self.shift -= 1;
        for exponent in 0..self.babies.len() as u32 {
            self.index(exponent);
        }
    }
}

/// The baby steps of a fixed base, answering the discrete logarithms of several targets.
///
/// [`discrete_log_shanks_steps`] grows its table while it searches, which is what a single
/// logarithm smaller than the order wants. A caller with several targets in the same group pays
/// instead for the whole table of about `sqrt(order)` baby steps once, and each query after it is
/// a walk of at most `sqrt(order)` giant steps over that table: Pohlig-Hellman queries it once per
/// digit of a prime power `pi**ri`, and the table is built once for the `ri` of them.
///
/// The memory is the one [`discrete_log_shanks_steps`] documents, so an order at or above
/// [`MAX_ORDER`] does not belong here either.
pub(crate) struct BabyStepTable<'r, R: ModRing> {
    /// The ring the base and the targets live in.
    ring: &'r R,
    /// `b**0, ..., b**(steps - 1)`, each mapped to its exponent.
    table: BabySteps<R::Elem>,
    /// Number of baby steps, which is also the stride of the giant steps.
    steps: u64,
    /// `b**-steps`: what a giant step multiplies the target by.
    stride: R::Elem,
    /// The order of the base, no exponent at or above it being a wanted answer.
    order: u64,
}

impl<'r, R: ModRing> BabyStepTable<'r, R> {
    /// The baby steps of `b`, whose order is `order`, `None` when `b` is not invertible.
    ///
    /// `order` must be below [`MAX_ORDER`], which is what bounds the memory of the table.
    pub(crate) fn new(ring: &'r R, b: &R::Elem, order: u64) -> Option<Self> {
        debug_assert!(order < MAX_ORDER, "an order the table cannot hold");
        let inverse = ring.invert(b)?;
        // An order of zero has a single exponent below it, as an order of one has.
        let order = order.max(1);
        // `steps**2 >= order`: the giant steps then cover every exponent below the order in at
        // most `steps` of them, as many as the baby steps cost.
        let steps = order.isqrt() + 1;

        let mut baby = ring.one();
        let mut table = BabySteps::new();
        table.push(baby.clone());
        for _ in 1..steps {
            ring.mul_assign(&mut baby, b);
            table.push(baby.clone());
        }

        Some(Self {
            ring,
            table,
            steps,
            stride: ring.pow(&inverse, &Integer::from(steps)),
            order,
        })
    }

    /// The smallest exponent `x` below the order with `b**x = a`, `None` when there is none.
    pub(crate) fn log(&self, a: &R::Elem) -> Option<u64> {
        // `giant = a * b**-offset`, so a baby step `b**j` equal to it means `a = b**(offset + j)`.
        let mut giant = a.clone();
        let mut offset = 0u64;
        loop {
            // This lookup covers the exponents `offset..offset + steps`, and the ones below
            // `offset` were covered by the previous ones: the smallest solution is the one found.
            if let Some(j) = self.table.get(&giant) {
                // A solution at or above the order only happens for an order that is not a
                // multiple of the real one, the last lookup reaching a little past it: the
                // smallest non-negative solution is the documented answer.
                return Some((offset + u64::from(j)) % self.order);
            }

            // Every exponent below the order has been covered.
            if offset + self.steps >= self.order {
                return None;
            }

            offset += self.steps;
            self.ring.mul_assign(&mut giant, &self.stride);
        }
    }
}

/// Baby-step giant-step algorithm for computing the discrete logarithm of `a` in base `b` modulo `n` (smallest non-negative integer `x` where `b**x = a (mod n)`).
///
/// The algorithm is a time-memory trade-off of the method of exhaustive search. It uses `O(sqrt(m))` memory, where `m` is the group order.
///
/// The baby steps and the giant steps are grown together, the stride of the giant steps following
/// the number of baby steps taken so far, and each giant step looks for the powers of the base met
/// until then: a logarithm `x` is found in about `sqrt(x)` steps with a table of about `sqrt(x)`
/// entries, instead of `sqrt(m)` of both whatever the logarithm is.
///
/// The baby steps stop at `sqrt(m) + 2` of them, and cost the residue each (8 bytes for a modulus
/// below `2**64`, 16 bytes and a GMP allocation above it) plus the 4 bytes of the slot indexing it,
/// of which two to four are held per baby step: 16 to 24 bytes each for a word-size modulus, so 16
/// to 24 MB for the largest order accepted. That figure is for word-size elements only. A modulus
/// above `2**64` holds each key in a GMP integer with an allocation of its own, which it does not
/// cover: the measured peak at the same order is 54 MB.
///
/// Orders of 10^12 or more are refused with [`Error::OutOfReach`], which says nothing about the
/// logarithm: it may well exist, and Pollard's rho will find it. The cap is a deliberately
/// conservative choice and not the size a table stops fitting in memory at, which those tens of
/// megabytes plainly do: it keeps the memory of a single call to a small fraction of what a process
/// can be expected to have, and the number of baby steps and of giant steps far inside a `u64`. Use
/// [`discrete_log_pollard_rho`](crate::discrete_log_pollard_rho), which builds no table at all, for
/// the larger orders.
///
/// A base that is not invertible modulo `n` is refused with [`Error::NotRelativelyPrime`], and a
/// modulus that is not positive with [`Error::InvalidModulus`]. Modulo 1 every residue is 0, so the
/// logarithm is 0.
///
/// If the order of the group is known, it can be passed as `order` to speed up the computation. An
/// `order` that is a proper multiple of the order of `b` still gives the smallest exponent here,
/// the exponents being covered in increasing order; that does not hold for every algorithm of this
/// crate, and [`crate::discrete_log_with_order`] documents what does.
pub fn discrete_log_shanks_steps(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    check_modulus(n)?;
    // Modulo 1 every residue is 0, so every exponent is a logarithm and 0 is the smallest one.
    if *n == 1 {
        return Ok(Integer::new());
    }
    let b = b.clone().modulo(n);
    let order = match order {
        Some(order) => order.clone(),
        None => n_order(&b, n)?,
    };

    if order >= MAX_ORDER {
        return Err(Error::OutOfReach);
    }
    // Below `MAX_ORDER` and not negative: exponents, the number of steps and the size of the
    // table all fit in a `u64`. Nothing is below a negative order.
    let Some(order) = order.to_u64() else {
        return Err(Error::LogDoesNotExist);
    };

    // The word-size ring whenever the modulus fits in one, GMP otherwise.
    match WordRing::from_modulus(n) {
        Some(ring) => {
            let (a, b) = (ring.from_integer(a), ring.from_integer(&b));
            solve(&ring, &a, &b, order)
        }
        None => {
            let ring = BigRing::new(n).expect("a positive modulus");
            let (a, b) = (ring.from_integer(a), ring.from_integer(&b));
            solve(&ring, &a, &b, order)
        }
    }
}

/// [`discrete_log_shanks_steps`] in `ring`, `a` and `b` being elements of it and `order` below
/// [`MAX_ORDER`].
fn solve<R: ModRing>(ring: &R, a: &R::Elem, b: &R::Elem, order: u64) -> Result<Integer, Error> {
    let inverse = ring.invert(b).ok_or(Error::NotRelativelyPrime)?;
    // An order of zero has a single exponent below it, as an order of one has.
    let order = order.max(1);
    // Baby steps stop at the square root of the order: past it a giant step covers more exponents
    // per multiplication than one more baby step would add.
    let last_baby = order.isqrt() + 1;

    let one = ring.one();
    let mut table = BabySteps::new();
    table.push(one.clone());

    // `b**top` is the last baby step in the table, and the table covers the exponents `0..=top`.
    let mut baby = one;
    let mut top = 0u64;
    // `giant = a * b**-offset`, so a baby step `b**j` equal to it means `a = b**(offset + j)`.
    let mut giant = a.clone();
    let mut offset = 0u64;
    // `b**-(top + 1)`, the stride of the giant steps: the exponents a giant step jumps over are
    // exactly the ones the table covers.
    let mut stride = inverse.clone();
    // What a full round of baby steps moves the stride by.
    let round_inverse = ring.pow(&inverse, &Integer::from(BABY_STEPS_PER_ROUND));

    loop {
        // This lookup covers the exponents `offset..=offset + top`, and the ones below `offset`
        // were covered by the previous ones: the smallest solution is the one found.
        if let Some(j) = table.get(&giant) {
            let x = offset + u64::from(j);
            // A solution above the order only happens for an order that is not a multiple of the
            // real one, the last lookup reaching a little past it: the documented result is the
            // smallest non-negative one.
            return Ok(Integer::from(x % order));
        }

        // Every exponent below the order has been covered.
        if offset + top + 1 >= order {
            return Err(Error::LogDoesNotExist);
        }

        // A giant step, over exactly the exponents the table covers.
        offset += top + 1;
        ring.mul_assign(&mut giant, &stride);

        // Baby steps for the next lookup, the stride following them.
        let mut added = 0;
        while added < BABY_STEPS_PER_ROUND && top < last_baby {
            ring.mul_assign(&mut baby, b);
            top += 1;
            table.push(baby.clone());
            added += 1;
        }
        if added == BABY_STEPS_PER_ROUND {
            ring.mul_assign(&mut stride, &round_inverse);
        } else {
            // The last round before the baby steps stop, and every round after it.
            for _ in 0..added {
                ring.mul_assign(&mut stride, &inverse);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rug::ops::Pow;

    use super::*;

    #[test]
    fn shanks_steps() {
        assert_eq!(
            discrete_log_shanks_steps(&442879.into(), &(Integer::from(7).pow(2)), &7.into(), None)
                .unwrap(),
            2
        );
        assert_eq!(
            discrete_log_shanks_steps(&874323.into(), &(Integer::from(5).pow(19)), &5.into(), None)
                .unwrap(),
            19
        );
        assert_eq!(
            discrete_log_shanks_steps(
                &6876342.into(),
                &(Integer::from(7).pow(71)),
                &7.into(),
                None
            )
            .unwrap(),
            71
        );
        assert_eq!(
            discrete_log_shanks_steps(
                &2456747.into(),
                &(Integer::from(3).pow(321)),
                &3.into(),
                None
            )
            .unwrap(),
            321
        );
    }

    #[test]
    fn order_out_of_reach() {
        // An order at the cap, and one far above it: the table is not built, and the refusal says
        // nothing about the logarithm, which here exists.
        for order in [
            Integer::from(MAX_ORDER),
            Integer::from(MAX_ORDER) * 1000u32,
            Integer::from(1) << 200u32,
        ] {
            assert_eq!(
                discrete_log_shanks_steps(
                    &Integer::from(265390227570863u64),
                    &Integer::from(184500076053622u64),
                    &2.into(),
                    Some(&order)
                ),
                Err(Error::OutOfReach),
                "order {order}"
            );
        }
        // One below it is accepted: the baby steps are interleaved with the giant steps, so a small
        // logarithm in a group of that size costs a handful of multiplications and no table to
        // speak of.
        let n = Integer::from(500_000_000_000u64).next_prime();
        let a = Integer::from(2).pow(11u32);
        assert_eq!(
            discrete_log_shanks_steps(&n, &a, &2.into(), Some(&Integer::from(MAX_ORDER - 1)))
                .unwrap(),
            11
        );
    }

    #[test]
    fn no_logarithm() {
        // The base is not invertible modulo `n`.
        assert_eq!(
            discrete_log_shanks_steps(&32942478.into(), &6.into(), &6.into(), Some(&1000.into())),
            Err(Error::NotRelativelyPrime)
        );
        // No power of the base is `a`.
        assert_eq!(
            discrete_log_shanks_steps(&442879.into(), &3.into(), &7.into(), None),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn smallest_solution() {
        // An order that is a multiple of the real one: every logarithm has several solutions
        // below it, and the smallest one is the answer. The repeated powers of the base also fill
        // the table with keys it already holds.
        for order in [3, 6, 9, 30] {
            for (a, x) in [(1, 0), (2, 1), (4, 2)] {
                assert_eq!(
                    discrete_log_shanks_steps(&7.into(), &a.into(), &2.into(), Some(&order.into()))
                        .unwrap(),
                    x,
                    "log of {a} in base 2 modulo 7, order {order}"
                );
            }
        }
    }

    #[test]
    fn small_logarithm_in_a_large_group() {
        // The baby steps of the whole order would be 700_000 entries here, where interleaving
        // stops a few steps in.
        let n = Integer::from(500_000_000_000u64).next_prime();
        let order = Integer::from(&n - 1u32);
        for x in [0u32, 1, 5, 1000] {
            let a = Integer::from(2).pow_mod(&Integer::from(x), &n).unwrap();
            assert_eq!(
                discrete_log_shanks_steps(&n, &a, &2.into(), Some(&order)).unwrap(),
                x
            );
        }
    }

    #[test]
    fn modulus_above_a_word() {
        // The moduli GMP is kept for, the order given since the real one is out of reach.
        let n = (Integer::from(1) << 70u32).next_prime();
        let a = Integer::from(3).pow_mod(&Integer::from(100), &n).unwrap();
        assert_eq!(
            discrete_log_shanks_steps(&n, &a, &3.into(), Some(&1000.into())).unwrap(),
            100
        );
        assert_eq!(
            discrete_log_shanks_steps(&n, &5.into(), &3.into(), Some(&1000.into())),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn invalid_modulus() {
        for n in [0, -1, -4] {
            assert_eq!(
                discrete_log_shanks_steps(&n.into(), &1.into(), &3.into(), None),
                Err(Error::InvalidModulus)
            );
        }
    }

    #[test]
    fn modulus_of_one() {
        // Modulo 1 every residue is 0, so the logarithm is 0 whatever `a`, `b` and the order are.
        for order in [None, Some(Integer::from(1)), Some(Integer::from(-5))] {
            assert_eq!(
                discrete_log_shanks_steps(&1.into(), &0.into(), &0.into(), order.as_ref()).unwrap(),
                0
            );
        }
        assert_eq!(
            discrete_log_shanks_steps(&1.into(), &5.into(), &7.into(), Some(&1000.into())).unwrap(),
            0
        );
    }

    #[test]
    fn shared_baby_steps_answer_several_queries() {
        // One table of the powers of 3 modulo 1000003, queried once per target: the same answers
        // as the interleaved search, which builds a table of its own for each of them.
        let n = Integer::from(1_000_003);
        let order = n_order(&3.into(), &n).unwrap();
        let ring = WordRing::from_modulus(&n).unwrap();
        let b = ring.from_integer(&3.into());
        let table = BabyStepTable::new(&ring, &b, order.to_u64().unwrap()).unwrap();
        for a in (1..1000u32).step_by(7) {
            let a = Integer::from(a);
            assert_eq!(
                table
                    .log(&ring.from_integer(&a))
                    .map(Integer::from)
                    .ok_or(Error::LogDoesNotExist),
                discrete_log_shanks_steps(&n, &a, &3.into(), Some(&order)),
                "log of {a} in base 3 modulo {n}"
            );
        }

        // An order that is a multiple of the real one: every logarithm has several solutions below
        // it, and the smallest one is the answer here too.
        let multiple = Integer::from(&order * 3u32);
        let table = BabyStepTable::new(&ring, &b, multiple.to_u64().unwrap()).unwrap();
        for a in (1..200u32).step_by(11) {
            let a = Integer::from(a);
            assert_eq!(
                table
                    .log(&ring.from_integer(&a))
                    .map(Integer::from)
                    .ok_or(Error::LogDoesNotExist),
                discrete_log_shanks_steps(&n, &a, &3.into(), Some(&multiple)),
                "log of {a} in base 3 modulo {n}, order {multiple}"
            );
        }

        // A base that is not invertible: no table at all, as the interleaved search refuses it.
        let ring = WordRing::from_modulus(&32942478.into()).unwrap();
        assert!(BabyStepTable::new(&ring, &ring.from_integer(&6.into()), 1000).is_none());
    }

    #[test]
    fn shared_baby_steps_above_a_word() {
        // The moduli GMP is kept for, the order given since the real one is out of reach: the keys
        // of the table are bignums, and the exponents covered still stop at the order.
        let n = (Integer::from(1) << 70u32).next_prime();
        let ring = BigRing::new(&n).unwrap();
        let b = ring.from_integer(&3.into());
        let table = BabyStepTable::new(&ring, &b, 1000).unwrap();
        for x in [0u64, 1, 2, 100, 999] {
            assert_eq!(table.log(&ring.pow(&b, &Integer::from(x))), Some(x));
        }
        assert_eq!(table.log(&ring.from_integer(&5.into())), None);
    }

    #[test]
    fn baby_steps_table() {
        let mut table = BabySteps::new();
        // Keys landing in the same slot of the smallest table: probing is what finds them again.
        let slot = |key: u64| hash_key(&key) >> (u64::BITS - INITIAL_SLOTS.trailing_zeros());
        let colliding: Vec<u64> = (1..100_000u64)
            .filter(|&k| slot(k) == slot(1))
            .take(3)
            .collect();
        assert_eq!(colliding.len(), 3);
        for &key in &colliding {
            table.push(key);
        }
        for (exponent, &key) in colliding.iter().enumerate() {
            assert_eq!(table.get(&key), Some(exponent as u32));
        }
        // A key already in the table keeps its exponent: the smallest one is the wanted one.
        table.push(colliding[0]);
        assert_eq!(table.get(&colliding[0]), Some(0));
        assert_eq!(table.get(&0), None);
    }

    #[test]
    fn baby_steps_table_growth() {
        // Every entry survives the rehashing of a table that doubled several times, keys of a size
        // GMP allocates for included.
        let mut table = BabySteps::new();
        let entries = 1000u32;
        for exponent in 0..entries {
            table.push(Integer::from(exponent) << 40u32);
        }
        assert!(table.slots.len() >= 2 * entries as usize);
        for exponent in 0..entries {
            assert_eq!(
                table.get(&(Integer::from(exponent) << 40u32)),
                Some(exponent)
            );
        }
        assert_eq!(table.get(&(Integer::from(entries) << 40u32)), None);
    }
}
