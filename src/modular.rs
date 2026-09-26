//! Modular arithmetic for moduli that fit in a machine word.
//!
//! GMP is built for numbers of any size: every operation is a call that checks sizes and loops
//! over limbs, and a division is expensive even on a single one. Almost every modulus met in
//! practice is below `2**64`, where a modular multiplication is one 128-bit multiplication and a
//! couple of shifts: [`WordRing`] keeps such a modulus in a `u64` and reduces with Montgomery
//! multiplication, [`BigRing`] falls back to GMP for the moduli that do not fit.
//!
//! Both implement [`ModRing`], so an algorithm written once over that trait runs natively on the
//! small moduli and on GMP on the large ones: `WordRing::from_modulus` returns the word-size ring
//! whenever the modulus fits in one, and `BigRing::new` is what the rest fall back to. The
//! `generic_over_any_modulus` test below is that pattern, whole.

use std::{cmp::Ordering, hash::Hash};

use rug::Integer;

/// Arithmetic in the ring of integers modulo a fixed positive modulus.
///
/// How an element is represented is up to the implementation ([`WordRing`] keeps its elements in
/// Montgomery form), so an element only means anything to the ring that produced it. Elements
/// are built with [`from_integer`](ModRing::from_integer), [`zero`](ModRing::zero) and
/// [`one`](ModRing::one), and read back with [`to_integer`](ModRing::to_integer).
///
/// Modulo 1 every element is 0, [`one`](ModRing::one) included.
// `from_integer` takes a receiver because the modulus belongs to the ring, not to the value: it
// converts into the ring, it does not build a ring.
#[allow(clippy::wrong_self_convention)]
pub trait ModRing {
    /// How an element of the ring is represented.
    ///
    /// The representation of a residue is unique: two elements are equal exactly when they are
    /// the same residue, and an element can be used as a hash key (what baby-step giant-step
    /// needs of it).
    type Elem: Clone + Eq + Hash;

    /// The element 0.
    fn zero(&self) -> Self::Elem;

    /// The element 1.
    fn one(&self) -> Self::Elem;

    /// The element `x mod n`, `x` being of any sign.
    fn from_integer(&self, x: &Integer) -> Self::Elem;

    /// The representative of `x` in `[0, n)`.
    fn to_integer(&self, x: &Self::Elem) -> Integer;

    /// The low word of the representative of `x` in `[0, n)`.
    ///
    /// An algorithm that branches on the value of an element, and not on however the ring happens
    /// to represent it (a random walk partitioning the group, for instance), reads it with this
    /// instead of [`to_integer`](ModRing::to_integer): no allocation, and no modulus of any size
    /// to divide by.
    fn residue_word(&self, x: &Self::Elem) -> u64 {
        self.to_integer(x).to_u64_wrapping()
    }

    /// `a + b`.
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    /// `a - b`.
    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    /// `a * b`.
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem;

    /// `a * a`.
    fn square(&self, a: &Self::Elem) -> Self::Elem;

    /// `a += b`.
    ///
    /// A loop that overwrites its accumulator uses the mutating operations instead of the ones
    /// above: a ring whose elements live on the heap ([`BigRing`]) reuses the buffer of `a`, where
    /// the by-value form allocates the result and frees the old value at every step. For a ring
    /// whose elements are words there is nothing to reuse, and the default is the whole
    /// implementation.
    #[inline]
    fn add_assign(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.add(a, b);
    }

    /// `a -= b`.
    #[inline]
    fn sub_assign(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.sub(a, b);
    }

    /// `a *= b`.
    #[inline]
    fn mul_assign(&self, a: &mut Self::Elem, b: &Self::Elem) {
        *a = self.mul(a, b);
    }

    /// `a *= a`.
    #[inline]
    fn square_assign(&self, a: &mut Self::Elem) {
        *a = self.square(a);
    }

    /// `a = a ** exponent`.
    ///
    /// # Panics
    ///
    /// If `exponent` is negative. A negative exponent is not an inverse power here: the two rings
    /// would answer it differently (the square and multiply ladder below reads the limbs of the
    /// absolute value, GMP's own powering wants an invertible base), so neither answers it. The
    /// check is unconditional, and the same one in both rings, so that a debug build and a release
    /// build of the same caller behave the same way.
    fn pow_assign(&self, a: &mut Self::Elem, exponent: &Integer) {
        check_exponent(exponent);
        square_and_multiply(self, a, exponent);
    }

    /// `a ** exponent`.
    ///
    /// # Panics
    ///
    /// If `exponent` is negative, as [`pow_assign`](ModRing::pow_assign) does.
    fn pow(&self, a: &Self::Elem, exponent: &Integer) -> Self::Elem {
        let mut result = a.clone();
        self.pow_assign(&mut result, exponent);
        result
    }

    /// The inverse of `a`, `None` when `a` and the modulus are not relatively prime.
    fn invert(&self, a: &Self::Elem) -> Option<Self::Elem>;
}

/// Refuses a negative exponent, the one thing both rings check before they raise anything.
///
/// `cmp0` is the sign of the number, which GMP keeps in the header of it, where a comparison with
/// the `Integer` 0 would be a call into GMP. The panic is out of line and cold, so what is left on
/// the path of every powering is a test of one word and a branch that is never taken.
///
/// Even that much inside [`ModRing::pow_assign`] used to cost Pollard's rho 1.4% of its
/// instructions, and none of it in the powering: the two functions that draw the multipliers and
/// the starting point of a walk were inlined into the walk itself, and growing them moved the
/// codegen of its inner loop by about one instruction per step. They are `#[inline(never)]` for that
/// reason, ten calls per search against millions of steps, and the check is then free.
#[inline]
fn check_exponent(exponent: &Integer) {
    if exponent.cmp0() == Ordering::Less {
        negative_exponent();
    }
}

/// The panic of [`check_exponent`], kept out of the code that checks.
#[cold]
#[inline(never)]
fn negative_exponent() -> ! {
    panic!("negative exponent");
}

/// `a = a ** exponent` by square and multiply, the exponent read limb by limb.
///
/// The limbs are those of the absolute value of the exponent, so a negative one would come out as
/// its opposite: [`ModRing::pow_assign`], the only way in, refuses it first.
///
/// Walking the exponent with [`Integer::get_bit`] is a call into GMP per bit, and each of them
/// costs more than the modular multiplication it decides: a candidate of the index calculus spent
/// half its instructions there. The limbs of the exponent are borrowed from GMP once instead, and
/// the loop shifts through them.
///
/// A fixed four-bit window was measured against this and lost by about 20% on every index calculus
/// instance: the exponents here are bounded by the group order, and over the 30 to 60 bits that
/// leaves, the fifteen multiplications building the table cost more than the ones it saves.
fn square_and_multiply<R: ModRing + ?Sized>(ring: &R, a: &mut R::Elem, exponent: &Integer) {
    // The limbs of the absolute value, least significant first, the most significant one not zero:
    // an exponent of 0 has none at all, and every element to it is 1.
    let Some((&high, rest)) = exponent.as_limbs().split_last() else {
        *a = ring.one();
        return;
    };
    // GMP chooses the width of its limbs when it is built, and it is not always the width of a
    // word here: the shifts below take it from the limb itself.
    let limb_bits = (std::mem::size_of_val(&high) * 8) as u32;

    // The leading 1 of the exponent needs neither the squarings above it nor a multiplication by
    // the base: the base itself is what the ladder starts from.
    let mut result = a.clone();
    for bit in (0..limb_bits - 1 - high.leading_zeros()).rev() {
        ring.square_assign(&mut result);
        if (high >> bit) & 1 != 0 {
            ring.mul_assign(&mut result, a);
        }
    }
    for &limb in rest.iter().rev() {
        for bit in (0..limb_bits).rev() {
            ring.square_assign(&mut result);
            if (limb >> bit) & 1 != 0 {
                ring.mul_assign(&mut result, a);
            }
        }
    }
    *a = result;
}

/// `a + b mod n`, `a` and `b` being below `n`.
#[inline]
fn add_mod(a: u64, b: u64, n: u64) -> u64 {
    let (sum, overflow) = a.overflowing_add(b);
    // The sum is below `2 * n`, so a single subtraction reduces it, and it is needed exactly when
    // the sum reached `n`: an overflow of the word means it did.
    if overflow || sum >= n {
        sum.wrapping_sub(n)
    } else {
        sum
    }
}

/// `a - b mod n`, `a` and `b` being below `n`.
#[inline]
fn sub_mod(a: u64, b: u64, n: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        // `a - b + n` is below `n`, the borrow and the addition cancel each other.
        a.wrapping_sub(b).wrapping_add(n)
    }
}

/// `n**-1 mod 2**64`, `n` being odd.
///
/// Newton's iteration doubles the number of correct bits at each step, starting from the three
/// bits an odd number is its own inverse on (`n**2 = 1 mod 8`).
fn inverse_mod_word(n: u64) -> u64 {
    debug_assert!(n % 2 == 1, "even modulus");
    let mut inverse = n;
    for _ in 0..5 {
        inverse = inverse.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inverse)));
    }
    debug_assert_eq!(n.wrapping_mul(inverse), 1);
    inverse
}

/// What Montgomery multiplication modulo an odd `n` needs, `R` being `2**64`.
#[derive(Clone, Debug)]
struct Montgomery {
    /// `-n**-1 mod R`: the factor that makes the low word of a reduction vanish.
    n_prime: u64,
    /// `R mod n`, the Montgomery form of 1.
    one: u64,
    /// `R**2 mod n`, what a plain residue is multiplied by to enter Montgomery form.
    r2: u64,
}

/// Arithmetic modulo a `u64`, the elements held in a single word.
///
/// An odd modulus is handled with Montgomery multiplication: elements are kept multiplied by
/// `2**64`, which turns a modular multiplication into two word multiplications instead of a
/// division. An even modulus has no such form (`2**64` is not invertible modulo it) and a plain
/// 128-bit remainder is used.
#[derive(Clone, Debug)]
pub struct WordRing {
    /// The modulus, never zero.
    modulus: u64,
    /// Montgomery constants, `None` when the modulus is even.
    montgomery: Option<Montgomery>,
}

// `from_u64` takes a receiver for the reason `ModRing::from_integer` does: the modulus belongs to
// the ring, not to the value.
#[allow(clippy::wrong_self_convention)]
impl WordRing {
    /// The ring of integers modulo `modulus`, `None` when the modulus is zero.
    pub fn new(modulus: u64) -> Option<Self> {
        if modulus == 0 {
            return None;
        }
        let montgomery = (modulus % 2 == 1).then(|| {
            // `R mod n`, without a 128-bit division.
            let one = (u64::MAX % modulus + 1) % modulus;
            // `R**2 mod n = (R mod n) * R`, reached by doubling once per bit of `R`.
            let mut r2 = one;
            for _ in 0..64 {
                r2 = add_mod(r2, r2, modulus);
            }
            Montgomery {
                n_prime: inverse_mod_word(modulus).wrapping_neg(),
                one,
                r2,
            }
        });
        Some(Self {
            modulus,
            montgomery,
        })
    }

    /// The ring of integers modulo `modulus`, `None` when the modulus is not positive or does not
    /// fit in a `u64`.
    ///
    /// This is the test an algorithm makes to choose between this ring and [`BigRing`].
    pub fn from_modulus(modulus: &Integer) -> Option<Self> {
        Self::new(modulus.to_u64()?)
    }

    /// The modulus.
    ///
    /// Nothing but the tests below reads it back: an algorithm passes the modulus to the ring and
    /// works with the elements it hands out ([`BigRing`] has no such accessor, nothing needing it).
    #[cfg(test)]
    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// The element `x mod n`.
    #[inline]
    pub fn from_u64(&self, x: u64) -> u64 {
        let x = x % self.modulus;
        match &self.montgomery {
            Some(montgomery) => self.redc(montgomery, u128::from(x) * u128::from(montgomery.r2)),
            None => x,
        }
    }

    /// The representative of `x` in `[0, n)`.
    #[inline]
    pub fn to_u64(&self, x: u64) -> u64 {
        match &self.montgomery {
            Some(montgomery) => self.redc(montgomery, u128::from(x)),
            None => x,
        }
    }

    /// `t / R mod n`, in `[0, n)`, `t` being below `n * R`.
    #[inline]
    fn redc(&self, montgomery: &Montgomery, t: u128) -> u64 {
        let low = t as u64;
        let high = (t >> 64) as u64;
        // `m * n` is what cancels the low word of `t`, making the division by `R` a shift.
        let m = low.wrapping_mul(montgomery.n_prime);
        let mn = u128::from(m) * u128::from(self.modulus);
        // Only the carry of the low words into the high ones is left to add: their sum is zero
        // modulo `R`, so it carries exactly when it is not zero at all.
        let (sum, carry) = high.overflowing_add((mn >> 64) as u64);
        let (sum, carry_again) = sum.overflowing_add(u64::from(low != 0));
        // `t + m * n` is below `2 * n * R`, so the quotient is below `2 * n`: one subtraction
        // reduces it, the carry out of the word counting as one `R` above the modulus.
        if carry || carry_again || sum >= self.modulus {
            sum.wrapping_sub(self.modulus)
        } else {
            sum
        }
    }
}

impl ModRing for WordRing {
    type Elem = u64;

    fn zero(&self) -> u64 {
        // The Montgomery form of 0 is 0.
        0
    }

    fn one(&self) -> u64 {
        match &self.montgomery {
            Some(montgomery) => montgomery.one,
            None => 1 % self.modulus,
        }
    }

    fn from_integer(&self, x: &Integer) -> u64 {
        // Entering the ring is not a hot path: an algorithm does it once per input.
        let residue = x.clone().modulo(&Integer::from(self.modulus));
        self.from_u64(
            residue
                .to_u64()
                .expect("a residue of a `u64` fits in a `u64`"),
        )
    }

    fn to_integer(&self, x: &u64) -> Integer {
        Integer::from(self.to_u64(*x))
    }

    #[inline]
    fn residue_word(&self, x: &u64) -> u64 {
        self.to_u64(*x)
    }

    #[inline]
    fn add(&self, a: &u64, b: &u64) -> u64 {
        add_mod(*a, *b, self.modulus)
    }

    #[inline]
    fn sub(&self, a: &u64, b: &u64) -> u64 {
        sub_mod(*a, *b, self.modulus)
    }

    #[inline]
    fn mul(&self, a: &u64, b: &u64) -> u64 {
        let product = u128::from(*a) * u128::from(*b);
        match &self.montgomery {
            Some(montgomery) => self.redc(montgomery, product),
            None => (product % u128::from(self.modulus)) as u64,
        }
    }

    #[inline]
    fn square(&self, a: &u64) -> u64 {
        self.mul(a, a)
    }

    fn invert(&self, a: &u64) -> Option<u64> {
        // An inversion is rare enough to be left to GMP's extended gcd.
        let inverse = Integer::from(self.to_u64(*a))
            .invert(&Integer::from(self.modulus))
            .ok()?;
        Some(
            self.from_u64(
                inverse
                    .to_u64()
                    .expect("a residue of a `u64` fits in a `u64`"),
            ),
        )
    }
}

/// Arithmetic modulo an [`Integer`], the fallback for the moduli [`WordRing`] cannot hold.
#[derive(Clone, Debug)]
pub struct BigRing {
    /// The modulus, always positive.
    modulus: Integer,
}

impl BigRing {
    /// The ring of integers modulo `modulus`, `None` when the modulus is not positive.
    pub fn new(modulus: &Integer) -> Option<Self> {
        (*modulus > 0).then(|| Self {
            modulus: modulus.clone(),
        })
    }
}

impl ModRing for BigRing {
    type Elem = Integer;

    fn zero(&self) -> Integer {
        Integer::new()
    }

    fn one(&self) -> Integer {
        Integer::from(1) % &self.modulus
    }

    fn from_integer(&self, x: &Integer) -> Integer {
        x.clone().modulo(&self.modulus)
    }

    fn to_integer(&self, x: &Integer) -> Integer {
        x.clone()
    }

    fn residue_word(&self, x: &Integer) -> u64 {
        // The low limb of the residue, without copying the rest of it.
        x.to_u64_wrapping()
    }

    fn add(&self, a: &Integer, b: &Integer) -> Integer {
        let mut sum = Integer::from(a + b);
        if sum >= self.modulus {
            sum -= &self.modulus;
        }
        sum
    }

    fn sub(&self, a: &Integer, b: &Integer) -> Integer {
        let mut difference = Integer::from(a - b);
        if difference < 0 {
            difference += &self.modulus;
        }
        difference
    }

    fn mul(&self, a: &Integer, b: &Integer) -> Integer {
        let mut product = Integer::from(a * b);
        product %= &self.modulus;
        product
    }

    fn square(&self, a: &Integer) -> Integer {
        let mut square = a.clone().square();
        square %= &self.modulus;
        square
    }

    // The mutating operations below are the reason this trait has them: GMP reduces into the
    // buffer of the left operand, so a loop of them allocates nothing, where every by-value
    // operation above allocates its result and frees what the caller replaces.

    fn add_assign(&self, a: &mut Integer, b: &Integer) {
        *a += b;
        if *a >= self.modulus {
            *a -= &self.modulus;
        }
    }

    fn sub_assign(&self, a: &mut Integer, b: &Integer) {
        *a -= b;
        if *a < 0 {
            *a += &self.modulus;
        }
    }

    fn mul_assign(&self, a: &mut Integer, b: &Integer) {
        *a *= b;
        *a %= &self.modulus;
    }

    fn square_assign(&self, a: &mut Integer) {
        a.square_mut();
        *a %= &self.modulus;
    }

    fn pow_assign(&self, a: &mut Integer, exponent: &Integer) {
        check_exponent(exponent);
        // GMP's own sliding window beats the default square and multiply.
        a.pow_mod_mut(exponent, &self.modulus)
            .expect("a non negative exponent");
    }

    fn invert(&self, a: &Integer) -> Option<Integer> {
        a.clone().invert(&self.modulus).ok()
    }
}

#[cfg(test)]
mod tests {
    use rug::rand::RandState;

    use super::*;

    #[test]
    #[should_panic(expected = "negative exponent")]
    fn word_ring_refuses_a_negative_exponent() {
        // The ladder reads the limbs of the absolute value, so it would answer `2**5` here: the
        // check is unconditional, so a release build refuses it as a debug build does.
        let ring = WordRing::new(587).unwrap();
        ring.pow(&ring.from_u64(2), &Integer::from(-5));
    }

    #[test]
    #[should_panic(expected = "negative exponent")]
    fn big_ring_refuses_a_negative_exponent() {
        // The same refusal, where GMP would otherwise panic on its own inside `pow_mod`.
        let ring = BigRing::new(&Integer::from(587)).unwrap();
        ring.pow(&ring.from_integer(&Integer::from(2)), &Integer::from(-5));
    }

    /// Moduli of every shape: 1, even, odd, powers of two, and the extremes of the word.
    const MODULI: [u64; 14] = [
        1,
        2,
        3,
        4,
        10,
        587,
        65536,
        4294967291,
        1 << 63,
        (1 << 63) + 1,
        u64::MAX - 2,
        u64::MAX - 1,
        u64::MAX,
        // The largest prime below `2**64`.
        18446744073709551557,
    ];

    /// Values to compute with, of every size up to a little above `2**64` and of both signs.
    fn values() -> Vec<Integer> {
        let mut rand = RandState::new();
        rand.seed(&Integer::from(0xD15Cu32));
        let mut values = vec![
            Integer::new(),
            Integer::from(1),
            Integer::from(-1),
            Integer::from(u64::MAX),
            Integer::from(u64::MAX) + 1u32,
            Integer::from(u64::MAX) * u64::MAX,
            -(Integer::from(u64::MAX) * 3u32),
        ];
        for _ in 0..8 {
            let x = Integer::from(Integer::random_bits(70, &mut rand));
            values.push(x.clone());
            values.push(-x);
        }
        values
    }

    /// Exponents of every length the ladder walks differently: none, a single bit, a few, the
    /// whole of one limb, and several limbs.
    fn exponents() -> Vec<Integer> {
        vec![
            Integer::new(),
            Integer::from(1),
            Integer::from(2),
            Integer::from(7),
            Integer::from(64),
            // The bit just below a word, then every bit of it.
            Integer::from(u64::MAX >> 1),
            Integer::from(u64::MAX),
            // One bit above a word, where a second limb starts, then two and three of them.
            Integer::from(u64::MAX) + 1u32,
            Integer::from(u64::MAX) * u64::MAX,
            (Integer::from(1) << 200u32) - 1u32,
        ]
    }

    /// Checks every operation of `ring` against the same computation done by GMP.
    fn check_ring<R: ModRing>(ring: &R, modulus: &Integer) {
        let values = values();
        assert_eq!(ring.to_integer(&ring.zero()), 0);
        assert_eq!(ring.to_integer(&ring.one()), Integer::from(1) % modulus);

        for a in &values {
            let ea = ring.from_integer(a);
            let ra = a.clone().modulo(modulus);
            assert_eq!(ring.to_integer(&ea), ra, "{a} mod {modulus}");
            assert_eq!(
                ring.to_integer(&ring.square(&ea)),
                ra.clone().square() % modulus,
                "{a}**2 mod {modulus}"
            );
            // The mutating operations compute the same thing as the ones returning a value: what
            // the hot loops rely on, and what the buffer reuse of `BigRing` could get wrong.
            let mut assigned = ea.clone();
            ring.square_assign(&mut assigned);
            assert_eq!(
                ring.to_integer(&assigned),
                ra.clone().square() % modulus,
                "{a}**2 mod {modulus} in place"
            );
            // Two elements are equal exactly when they are the same residue.
            assert!(ea == ring.from_integer(&ra));
            assert!(
                ea == ring.from_integer(&(ra.clone() + modulus)),
                "{a} + {modulus}"
            );

            for exponent in exponents() {
                let expected = ra.clone().pow_mod(&exponent, modulus).unwrap();
                assert_eq!(
                    ring.to_integer(&ring.pow(&ea, &exponent)),
                    expected,
                    "{a}**{exponent} mod {modulus}"
                );
                let mut assigned = ea.clone();
                ring.pow_assign(&mut assigned, &exponent);
                assert_eq!(
                    ring.to_integer(&assigned),
                    expected,
                    "{a}**{exponent} mod {modulus} in place"
                );
            }

            match ring.invert(&ea) {
                Some(inverse) => assert_eq!(
                    ring.to_integer(&inverse),
                    ra.clone().invert(modulus).unwrap(),
                    "1/{a} mod {modulus}"
                ),
                None => assert!(ra.clone().invert(modulus).is_err(), "1/{a} mod {modulus}"),
            }

            for b in &values {
                let eb = ring.from_integer(b);
                let rb = b.clone().modulo(modulus);
                for (name, expected, by_value, assign) in [
                    ("+", (ra.clone() + &rb) % modulus, ring.add(&ea, &eb), {
                        let mut assigned = ea.clone();
                        ring.add_assign(&mut assigned, &eb);
                        assigned
                    }),
                    (
                        "-",
                        (ra.clone() - &rb).modulo(modulus),
                        ring.sub(&ea, &eb),
                        {
                            let mut assigned = ea.clone();
                            ring.sub_assign(&mut assigned, &eb);
                            assigned
                        },
                    ),
                    ("*", (ra.clone() * &rb) % modulus, ring.mul(&ea, &eb), {
                        let mut assigned = ea.clone();
                        ring.mul_assign(&mut assigned, &eb);
                        assigned
                    }),
                ] {
                    assert_eq!(
                        ring.to_integer(&by_value),
                        expected,
                        "{a} {name} {b} mod {modulus}"
                    );
                    assert_eq!(
                        ring.to_integer(&assign),
                        expected,
                        "{a} {name} {b} mod {modulus} in place"
                    );
                }
            }
        }
    }

    #[test]
    fn word_ring() {
        for modulus in MODULI {
            let ring = WordRing::new(modulus).unwrap();
            assert_eq!(ring.modulus(), modulus);
            check_ring(&ring, &Integer::from(modulus));
        }
    }

    #[test]
    fn word_ring_words() {
        // `from_u64` and `to_u64` are the same conversions, without going through GMP.
        for modulus in MODULI {
            let ring = WordRing::new(modulus).unwrap();
            for x in [0, 1, 2, 587, u64::MAX / 2, u64::MAX - 1, u64::MAX] {
                assert_eq!(ring.to_u64(ring.from_u64(x)), x % modulus);
                assert_eq!(ring.from_u64(x), ring.from_integer(&Integer::from(x)));
            }
        }
    }

    #[test]
    fn residue_words() {
        // The same word for the same residue, whatever the ring: what a walk partitioning a group
        // by the value of its elements relies on.
        for modulus in MODULI {
            let ring = WordRing::new(modulus).unwrap();
            let big = BigRing::new(&Integer::from(modulus)).unwrap();
            for x in &values() {
                let residue = x.clone().modulo(&Integer::from(modulus));
                assert_eq!(
                    ring.residue_word(&ring.from_integer(x)),
                    residue.to_u64_wrapping(),
                    "{x} mod {modulus}"
                );
                assert_eq!(
                    big.residue_word(&big.from_integer(x)),
                    residue.to_u64_wrapping(),
                    "{x} mod {modulus}"
                );
            }
        }
        // A modulus no word can hold: only the low limb of the residue is kept.
        let modulus = Integer::from(u64::MAX) * u64::MAX + 1u32;
        let big = BigRing::new(&modulus).unwrap();
        let x = Integer::from(u64::MAX) + 1u32;
        assert_eq!(big.residue_word(&big.from_integer(&x)), 0);
    }

    #[test]
    fn big_ring() {
        for modulus in MODULI {
            let modulus = Integer::from(modulus);
            check_ring(&BigRing::new(&modulus).unwrap(), &modulus);
        }
        // A modulus no word can hold.
        let modulus = Integer::from(u64::MAX) * u64::MAX + 1u32;
        check_ring(&BigRing::new(&modulus).unwrap(), &modulus);
    }

    #[test]
    fn invalid_moduli() {
        assert!(WordRing::new(0).is_none());
        assert!(WordRing::from_modulus(&Integer::new()).is_none());
        assert!(WordRing::from_modulus(&Integer::from(-7)).is_none());
        // Just above the largest modulus a word holds.
        assert!(WordRing::from_modulus(&(Integer::from(u64::MAX) + 1u32)).is_none());
        assert!(WordRing::from_modulus(&Integer::from(u64::MAX)).is_some());

        assert!(BigRing::new(&Integer::new()).is_none());
        assert!(BigRing::new(&Integer::from(-7)).is_none());
    }

    #[test]
    fn generic_over_any_modulus() {
        // What an algorithm of the crate is written like: one body over the trait, and the ring
        // chosen from the size of the modulus.
        fn cube_of_3<R: ModRing>(ring: &R) -> Integer {
            let three = ring.from_integer(&Integer::from(3));
            ring.to_integer(&ring.mul(&three, &ring.square(&three)))
        }

        for modulus in [Integer::from(7u32), Integer::from(u64::MAX) * u64::MAX] {
            let cube = match WordRing::from_modulus(&modulus) {
                Some(ring) => cube_of_3(&ring),
                None => cube_of_3(&BigRing::new(&modulus).unwrap()),
            };
            assert_eq!(cube, Integer::from(27) % &modulus);
        }
    }

    #[test]
    fn large_exponents() {
        // Exponents larger than the modulus, where the square and multiply loop is long.
        let exponent = Integer::from(u64::MAX) * 3u32;
        for modulus in MODULI {
            let ring = WordRing::new(modulus).unwrap();
            let big = BigRing::new(&Integer::from(modulus)).unwrap();
            let base = Integer::from(12345);
            assert_eq!(
                ring.to_integer(&ring.pow(&ring.from_integer(&base), &exponent)),
                big.to_integer(&big.pow(&big.from_integer(&base), &exponent)),
                "12345**{exponent} mod {modulus}"
            );
        }
    }
}
